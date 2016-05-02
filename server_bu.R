
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com

library(shiny)

shinyServer(
  function(input, output, clientData, session) {
    
    #------------------------------------------------------
    #------------------------------------------------------
    # UPDATE THE UI DATA
    #------------------------------------------------------
    #------------------------------------------------------
    
    observe({
      
      rs <- Data()$reporter
      rs1 <- Data()$gene
      
      reps <- na.omit(unique(factor(rs$line)))
      gens <- na.omit(factor(unique(rs1$Gene_ID)))
      
      # Genotype list
      s_options <- list()
      for(r in reps) s_options[[r]] <- r
      updateSelectInput(session, "ref_reps", choices = s_options)  
      
      # Genes
      g_options <- list()
      for(g in gens) g_options[[g]] <- g
      updateSelectInput(session, "ref_genes", choices = g_options)  
      
      # Genotype check box
      updateSelectInput(session, "to_plot", choices = s_options, selected = 2)  
  
    })
    
    
    #------------------------------------------------------
    #------------------------------------------------------
    # LOAD A DIRECTORY
    #------------------------------------------------------
    #------------------------------------------------------    
    observeEvent(
      ignoreNULL = TRUE,
      eventExpr = {
        input$directory
      },
      handlerExpr = {
        if (input$directory > 0) {
          # condition prevents handler execution on initial app launch
          
          # launch the directory selection dialog with initial path read from the widget
          path = choose.dir(default = readDirectoryInput(session, 'directory'))
          
          # update the widget value
          updateDirectoryInput(session, 'directory', value = path)
        }
      }
    )
    
    
    #------------------------------------------------------
    #------------------------------------------------------
    # LOAD THE USER DATA
    #------------------------------------------------------
    #------------------------------------------------------
    
    Data <- reactive({
      
      # Load two datafiles
      inGene <- input$gene_file
      pathData <- readDirectoryInput(session, 'directory') # this is a specific function to load the directory
      
      data <- list()

      if (is.null(pathData) | is.null(inGene)) return(NULL)
    
      # Attach the gene informations
      data$gene <- read.table(inGene$datapath, header = T) 
      
      message(pathData)
      
      # Attach the reporter informations
      list.f <- list.files(pathData)
      rs <- NULL
      for(f in list.f){
        name <- gsub(".xlsx", "", f)
        for(i in 1:20){
          tryCatch({
            temp <- read_excel(paste0(pathData, f), sheet = i)
            temp <- temp[!is.na(temp[,1]),]
            rs <- rbind(rs, data.frame(line=name, root=i, cell_type=temp$Label, value=temp[["Average flourescence"]]))
          },warning = function(w) {
          }, error = function(e) {
          })
        }
      }
      remove(temp, i, f, list.f, name)
      
      # Average the data by line, root, cell type
      mean_data <- ddply(rs, .(line, root, cell_type), summarise, value=mean(value))
      mean_data <- mean_data[!is.na(mean_data$value),]
      
      # Reshape the data to have them in the proper form for the analysis
      reporter <- dcast(mean_data, line + root ~ cell_type)
      
      data$reporter <- reporter
      
      passed <<- FALSE
      
      return(data)
      
    })  
  
    
    #------------------------------------------------------
    #------------------------------------------------------
    # PROCESS THE DATA
    #------------------------------------------------------
    #------------------------------------------------------
  
     
    Results <- reactive({
      

      if(input$runROOTEXP == 0){return()}
      
      #rs <- read.csv("~/Desktop/root_data.csv", header = T) 
      rs <- Data()$reporter
      #rs <- rs[,-1]
      rs <- na.omit(rs)
      
      rs1 <- rs
      p.list <- as.character(unique(rs$line))
      tissues <- colnames(rs)[-c(1,2)]
      
      
# Normalize the dataset
      #rs <- normalize(rs)
      
      
      # Aggregate the data      
      rs.melt <<- melt(rs, id.vars =c("line", "root"))
      rs.agg <<- ddply(rs.melt, .(line, variable), summarize, value=mean(value))
      rs.agg$type <- F
      
      # Get the difference with the wild type
      rs.agg$wt <- rep(rs.agg$value[rs.agg$line == input$ref_reps], length(p.list))
      rs.agg$diff <- rs.agg$value - rs.agg$wt
      rs.agg$x <- as.numeric(rs.agg$variable)
      
# ----------------------------------------------------------------------------------
# ----- PAIRWISE ANOVA COMPARISONS AND MANOVA ANALYSIS   ---------------------------
# ----------------------------------------------------------------------------------            
      
    
      # if(!passed){
      isolate({
        passed <<- TRUE
        
        # Create a table that will contain the anova results
        message("Computing the differences between tissues")
        l1 <- length(p.list)
        l2 <- l1^2
        l3 <- (((l1 * l1) - l1) / 2) * length(tissues)
        aov.results <<- data.frame(genotype_1 = character(l3), genotype_2=character(l3), tissue = character(l3), pvalue = numeric(l3), stringsAsFactors = F)
        maov.results <<- matrix(0, ncol = l1, nrow = l1); colnames(maov.results) <- p.list; rownames(maov.results) <- p.list
        i <- 1
        k <- 1
        for(p in p.list){
          j <- 1
          for(p1 in p.list){
            if(p != p1 & maov.results[j,i] == 0){
              
                # MANOVA analysis to compare the plants
                temp <- rs[rs$line == p | rs$line == p1,]
                fit <- manova(cbind(lateralrootcap, columella, QC, epidermis, cortex, endodermis) ~ line,  temp) # Stele was removed
                maov.results[j,i] <- round(summary(fit, test = "Wilks")$stats[1,6], 5)
                maov.results[i,j] <- maov.results[j,i]
                
                # Reponse for each tissue
                for(ti in tissues){
                  temp <- rs.melt[rs.melt$line == p | rs.melt$line == p1,]
                  temp1 <- temp[temp$variable == ti,]
                  if(nrow(temp1) > 3){
                    fit <- aov(value ~ line, data=temp1)
                    aov.results[k,] <- c(p1, p, ti, round(summary(fit)[[1]][1,5], 5))
                  }
                  else{ # In case there is noit enough data to make the anova
                    aov.results[k,] <- 1
                  }
                  k <- k+1
                }
                
            }
            j <- j+1
          }
          i <- i+1
        }

      # }
        
      })
        
# ----------------------------------------------------------------------------------
# ----- LDA ANALYSIS ---------------------------------------------------------------
# ----------------------------------------------------------------------------------            
      
      # Make the LDA analysis on the reporter dataset.
      temp <- rs[,-2]
      
      fit <- lda(line ~ ., data=temp)
      fit.p <- predict(fit, newdata = temp[,-1])
      rs <- cbind(rs, fit.p$x)

      
# ----------------------------------------------------------------------------------
# ------ ROOT MAP CREATION ---------------------------------------------------------
# ----------------------------------------------------------------------------------
      
      
      message("Creating the root maps")
      
      # Open the SVG as an XML file

# Create a color palette for the gene expression
      temp1 <- rs.agg[rs.agg$line == input$ref_reps | rs.agg$line == input$to_plot,]
      co <- temp1$value
      co <- co - min(co)
      temp1$col <- round((co + 1)*10)
      my.col <- colorRampPalette(c("yellow", "red"))(diff(range(temp1$col))+1)
      
# Make the Root Map for the wild type

      # Modifyt the SVG for each tissue
      temp <- temp1[temp1$line == input$ref_reps ,]
      root1 <- xmlParse("./www/root.svg")
      i <- 1
      for(n in tissues){
        node = xpathApply(root1, paste0("//*[@id='",n,"']/*"))
        sapply(node, function(x) {
          oldstyle <- xmlAttrs(x)
          removeAttributes(x, "style")
          xmlAttrs(x)<-c(gsub("fill:#fff", paste0("fill:",my.col[temp$col[temp$variable == n]]),
                oldstyle['style']))    
        })
        i <- i+1
      }
      # Save the new file
      saveXML(root1, file='./www/root_1.svg')
      
      
# Make the Root Map for the mutant
      
      # Modifyt the SVG for each tissue
      root2 <- xmlParse("./www/root.svg")
      temp <- temp1[temp1$line == input$to_plot ,]
      for(n in tissues){
        node = xpathApply(root2, paste0("//*[@id='",n,"']/*"))
        sapply(node, function(x) {
          oldstyle <- xmlAttrs(x)
          removeAttributes(x, "style")
          xmlAttrs(x)<-c(gsub("fill:#fff", paste0("fill:",my.col[temp$col[temp$variable == n]]),
                              oldstyle['style']))    
        })
      }
      # Save the new file
      saveXML(root2, file='./www/root_2.svg')
      
      
# Make the Root Map for the mutant-wildtype
      
      # Create a color palette for the gene expression
      temp <- rs.agg[rs.agg$line == input$to_plot,]
      temp$sig <- aov.results$pvalue[(aov.results$genotype_1 == input$to_plot & aov.results$genotype_2 == input$ref_reps) | 
                     (aov.results$genotype_1 == input$ref_reps & aov.results$genotype_2 == input$to_plot)] < 0.05
      temp <- temp[temp$sig == T,]
      if(nrow(temp) > 0){
        co <- round(temp$diff)
        if(min(co) < 0) co <- co + abs(min(co))
        if(min(co) > 0) co <- co - abs(min(co))
        temp$col <- co + 1
        my.col <- colorRampPalette(c("blue", "red"))(diff(range(co))+1)
        
        ti <- unique(temp$variable)
  
        # Modifyt the SVG for each tissue
        root_diff <- xmlParse("./www/root.svg")
        for(n in ti){
          node = xpathApply(root_diff, paste0("//*[@id='",n,"']/*"))
          sapply(node, function(x) {
            oldstyle <- xmlAttrs(x)
            removeAttributes(x, "style")
            xmlAttrs(x)<-c(gsub("fill:#fff", paste0("fill:",my.col[temp$col[temp$variable == n]]),
                                oldstyle['style']))    
          })
        }
      }else{
        root_diff <- xmlParse("./www/root.svg")
      }
      # Save the new file
      saveXML(root_diff, file='./www/root_diff.svg')      
      

      # Use the command line to convert the svg to png for further use. 
      # This is the step that takes time
      root1.bm <- rsvg("./www/root_1.svg"); writePNG(root1.bm, "./www/root_1.png")
      root2.bm <- rsvg("./www/root_2.svg"); writePNG(root2.bm, "./www/root_2.png")
      rootdiff.bm <- rsvg("./www/root_diff.svg"); writePNG(rootdiff.bm, "./www/root_diff.png")
      
      
      
      # Return the results
      
      results <<- list()
      results$data <<- rs1
      results$norm <<- rs
      results$agg <<- rs.agg
      results$maov <<- maov.results
      results$aov <<- maov.results
      results$p.list <<- p.list
      results$tissues <<- tissues
      #results$mean.diff <<- mean.diff
      results$nameWt <<- input$ref_reps
      results$nameMt <<- input$to_plot
      results$root_1 <<- './www/root_1.png'
      results$root_2 <<- './www/root_2.png'
      results$root_diff <<- './www/root_diff.png'
      
      return(results)  
      
    })
    
    
    
  
    # Plot the different distributions
    histoPlot <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      rs <- Results()$data
      plot1 <- ggplot(rs, aes(x=value, fill=plant, colour=plant)) +
        geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
        theme_bw() + 
        ggtitle("Original data")
      
      rs <- Results()$norm
      plot2 <- ggplot(rs, aes(x=value, fill=plant, colour=plant)) +
        geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
        theme_bw() +
        ggtitle("Normalized data")
      
      grid.arrange(plot1, plot2, nrow=1, ncol=2)
    }
    
    output$histoPlot <- renderPlot({
      print(histoPlot())
    })
    
    output$downloadPlot <- downloadHandler(
      filename = "distribution.png",
      content = function(file) {
        png(file, width = input$plot_width, height=input$plot_height)
        histoPlot()
        dev.off()
      })  
    
    
    
    
    # Plot the different distributions
    diffPlot <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      rs.agg <- Results()$agg
      rs.agg$mean.diff <- Results()$mean.diff
      tissues <- Results()$tissues

      plot1 <- ggplot(data = rs.agg, aes(x, diff, colour=plant)) + 
        geom_rect(aes(xmin=min(x), xmax=max(x), 
                      ymin=-mean.diff, ymax=mean.diff), 
                  fill="lightgrey", color="lightgrey", alpha=0.5) +
        geom_point(data=rs.agg[rs.agg$type == 1,], 
                   aes(x, diff, colour=plant), size=6, pch=1) +
        geom_line(size=1) + 
        geom_point(size=3) +
        scale_x_discrete(breaks= c(1:length(tissues)),labels=tissues) +
        theme_bw()
      plot1

    }
    
    output$diffPlot <- renderPlot({
      if(input$runROOTEXP == 0){return()}
      
      print(diffPlot())
    })
    
    output$downloadPlot <- downloadHandler(
      filename = "difference.png",
      content = function(file) {
        png(file, width = input$plot_width, height=input$plot_height)
        diffPlot()
        dev.off()
      })     
    
    
    
    
    # Plot the different distributions
    diffPlotSingle <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      rs.agg <- Results()$agg
      rs.agg$mean.diff <- Results()$mean.diff
      tissues <- Results()$tissues
      
      plot1 <- ggplot(rs.agg, aes(x, diff, colour=plant)) + 
        geom_rect(aes(xmin=min(rs.agg$x), xmax=max(rs.agg$x), 
                      ymin=-mean.diff, ymax=mean.diff), 
                  fill="lightgrey", color="lightgrey", alpha=0.5) +
        geom_line(size=0.5) + 
        geom_line(data = rs.agg[rs.agg$plant == input$to_plot,], aes(x, diff), colour="#FF7C00", size=1) + 
        geom_point(data=rs.agg[rs.agg$type == 1 & rs.agg$plant == input$to_plot,], aes(x, diff), colour="#FF7C00",size=6, pch=1) +
        geom_point(size=2) +
        geom_point(data=rs.agg[rs.agg$plant == input$to_plot,], aes(x, diff), colour="#FF7C00", size=3) +
        scale_x_discrete(breaks= c(1:length(tissues)),labels=tissues) +
        scale_colour_grey()+
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
      
      plot1
      
    }
    output$diffPlotSingle <- renderPlot({
      if(input$runROOTEXP == 0){return()}
      
      print(diffPlotSingle())
    })    
    
    
    # Plot the different distributions
    heatmap <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      temp <- Results()$maov
      temp[temp > 0.05] <- 0.06
      temp[temp == 0] <- 0.06
      ramp <- colorRamp(c("white", "lightblue",  "red"))
      plot1 <- levelplot(temp,  scales=list(x=list(rot=45)), col.regions=rgb(ramp(seq(1, 0, length = 100)), max = 255), xlab = "", ylab = "")
      plot1
      
      temp2 <- Results()$norm
      plot2 <- ggplot(temp2, aes(LD1, LD2, colour=line)) + 
        geom_point(size=4) +
        theme_bw() + 
        stat_ellipse()
      
      
    }
    output$heatmap <- renderPlot({
      print(heatmap())
    })   
 
    
    
    
    
    output$root_1 <- renderImage({
      
      if (input$runROOTEXP == 0) { return(" Please run analysis ")}
      
      filename <- normalizePath(file.path(Results()$root_1))
      
      # Return a list containing the filename and alt text
      list(src = filename, width=150)     
      
    }, deleteFile = FALSE
    )
    
    
    output$root_2 <- renderImage({
      
      if (input$runROOTEXP == 0) { return(" Please run analysis ")}

      filename <- normalizePath(file.path(Results()$root_2))
      message(filename)
      # Return a list containing the filename and alt text
      list(src = filename, width=150)     
      
    }, deleteFile = FALSE
    )

    
    output$root_diff <- renderImage({
      
      if (input$runROOTEXP == 0) { return(" Please run analysis ")}
      
      filename <- normalizePath(file.path(Results()$root_diff))
      
      # Return a list containing the filename and alt text
      list(src = filename, width=150)     
      
    }, deleteFile = FALSE
    )  
    
    output$nameWt <- renderText({ 
      if(input$runROOTEXP == 0){return(" Please run analysis ")}
      
      Results()$nameWt
    })    
    
    output$nameMt <- renderText({ 
      if(input$runROOTEXP == 0){return(" Please run analysis ")}
      Results()$nameMt
    })    
    
    output$nameDiff <- renderText({ 
      if(input$runROOTEXP == 0){return(" Please run analysis ")}
      paste0(Results()$nameMt, " - ", Results()$nameWt)
    })    
    
    

})
