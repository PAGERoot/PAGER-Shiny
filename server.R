
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com

library(shiny)

gene <<- NULL


shinyServer(
  function(input, output, clientData, session) {

    #------------------------------------------------------
    # CHOOSE DIRECTORY
    #------------------------------------------------------
#     volumes <- getVolumes() #c('R Installation'=R.home())
#     shinyDirChoose(input, 'directory', roots=volumes, session=session)
# 
#     observe({
#       dir_label <- parseDirPath(volumes, input$directory)
#       updateTextInput(session, "dirPath", value = dir_label)
#     })
    
        
    observe({
            
      if(input$load_data == 0){return()}

      withProgress(message = 'WORKING:', value = 0, {
      
          

        #------------------------------------------------------
        # LOAD THE USER DATA
        #------------------------------------------------------

        incProgress(1/4, detail = "Loading the data")
        

        if(input$use_example){ # Use the exmple dataset
          gene <<- read.table("www/GeneExpression.txt", header = T)   
          rs <<- read.csv("www/reporter-data.csv", header = T) 
          message(str(gene))
          message("-----------------")
        }else{
          # Load two datafiles
          inGene <- input$gene_file
          #readDirectoryInput(session, 'directory') # this is a specific function to load the directory
          # pathData <- paste0(input$dirPath, "/")
          inReporter <- input$rep_file
          
          if (is.null(inReporter)) return(NULL)
          if(!is.null(inGene)){
            gene <<- read.table(inGene$datapath, header = T)   
          }
          # Attach the gene informations
          
  #         # Attach the reporter informations
  #         list.f <- list.files(pathData)
  #         rs <- NULL
  #         for(f in list.f){
  #           name <- gsub(".xlsx", "", f)
  #           for(i in 1:20){
  #             tryCatch({
  #               temp <- read_excel(paste0(pathData, f), sheet = i)
  #               temp <- temp[!is.na(temp[,1]),]
  #               # Normalize the fluorescence data
  #               fluo <- scale(temp[["Average flourescence"]])
  #               rs <- rbind(rs, data.frame(line=name, root=i, cell_type=temp$Label, value=fluo))
  #             },warning = function(w) {
  #             }, error = function(e) {
  #             })
  #           }
  #         }
  #         remove(temp, i, f, list.f, name)
  #                 
          rs <<- read.csv(inReporter$datapath, header = T)   
        }
        # Average the data by line, root, cell type
        mean_data <- ddply(rs, .(line, root, cell_type), summarise, value=mean(value))
        mean_data <- mean_data[!is.na(mean_data$value),]
        
        # Reshape the data to have them in the proper form for the analysis
        reporter <- dcast(mean_data, line + root ~ cell_type)
        
        reporter <<- reporter      
      
      
      
      
      #------------------------------------------------------
      # PROCESS THE DATA
      #------------------------------------------------------

      
      # rs <- read.csv("~/Desktop/root_data.csv", header = T) 
      rs <- reporter
      rs <- na.omit(rs)
      
      p.list <- as.character(unique(rs$line))
      tissues <- colnames(rs)[-c(1,2)]
      
      
      # Aggregate the data      
      rs.melt <- melt(rs, id.vars =c("line", "root"))
      
      
      # ----------------------------------------------------------------------------------
      # ----- PAIRWISE ANOVA COMPARISONS AND MANOVA ANALYSIS   ---------------------------
      # ----------------------------------------------------------------------------------            
      
      incProgress(1/4, detail = "Performing the MANOVA")
      
      # Create a table that will contain the anova results
      l1 <- length(p.list)
      l2 <- l1^2
      l3 <- (((l1 * l1) - l1) / 2) * length(tissues)
      
      aov.results <<- data.frame(genotype_1 = character(l3), 
                                 genotype_2=character(l3), 
                                 tissue = character(l3), 
                                 pvalue = numeric(l3), 
                                 stringsAsFactors = F)
      
      maov.results <<- matrix(0, ncol = l1, nrow = l1) 
      colnames(maov.results) <- p.list 
      rownames(maov.results) <- p.list
      
      i <- 1
      k <- 1
      for(p in p.list){
        j <- 1
        for(p1 in p.list){
          if(p != p1 & maov.results[j,i] == 0){
            
            #print(paste0(p, " / ", p1))
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
      
      # -------------------------------------------------------
      # ----- LDA ANALYSIS ------------------------------------
      # -------------------------------------------------------            
      
      incProgress(1/4, detail = "Performing the LDA")
      
      # Make the LDA analysis on the reporter dataset.
      temp <- rs[,-2]
      fit <- lda(line ~ ., data=temp, CV=F)
      
      # Get the accuracy of the prediction
      fit.p <- predict(fit, newdata = temp[,-1])
      ct <- table(temp$line, fit.p$class)
      diag(prop.table(ct, 1))
      sum(diag(prop.table(ct)))
      
      #gene <- read.csv("~/Desktop/GeneExpression.txt", sep="\t")
      message(str(gene))
      if(!is.null(gene)){
        message("---- Prediction")
        prediction <- predict(fit, newdata=gene[,-1])$class
        gene <<- cbind(gene, prediction)
      }
      
      
      rs <- cbind(rs, fit.p$x)
      

      lda.fit <<- fit
      rs <<- rs
      rs.melt <<- rs.melt
      tissues <<- tissues
      rs.aov <<- aov.results
      rs.maov <<- maov.results
      p.list <<- p.list
      rs.agg.short <<- ddply(rs.melt, .(line, variable), summarize, value=mean(value))
      
      
      #------------------------------------------------------
      # ----- UPDATE THE UI DATA ----------------------------
      #------------------------------------------------------

      incProgress(1/4, detail = "Updating the UI")
      
      reps <- na.omit(unique(factor(reporter$line)))
      reps <- na.omit(unique(factor(rs$line)))
      
      # Genotype list
      s_options <- list()
      for(r in reps) s_options[[r]] <- r
      updateSelectInput(session, "ref_reps", choices = s_options)  
      
      # Genotype check box
      updateSelectInput(session, "to_plot", choices = s_options, selected=s_options[2])  
      
      # Genes
      if(!is.null(gene)) {
        gens <- na.omit(factor(unique(gene$Gene_ID)))
        g_options <- list()
        for(g in gens) g_options[[g]] <- g
        updateSelectInput(session, "ref_genes", choices = g_options)        
      }
      

      })
      
      
      
    })
    
    
    
  # ----------------------------------------------------------------------------------
  # ------ UPDATE UI TO PREVENT IDENTIQUE SELECTIONS           ------------------
  # ----------------------------------------------------------------------------------
    
    observe({
      
      if(input$runROOTEXP == 0){return()}
      
      sel <- input$to_plot
      reps <- na.omit(unique(factor(reporter$line)))

      reps1 <- reps[-match(input$ref_reps, reps)]
      
      # Genotype check box
      s_options <- list()
      for(r in reps1) s_options[[r]] <- r
      updateSelectInput(session, "to_plot", choices = s_options, selected=sel)  
    })
    
    observe({
      
      if(input$runROOTEXP == 0){return()}
      
      sel <- input$ref_reps
      reps <- na.omit(unique(factor(reporter$line)))

      reps2 <- reps[-match(input$to_plot, reps)]
      
      # Genotype list
      s_options <- list()
      for(r in reps2) s_options[[r]] <- r
      updateSelectInput(session, "ref_reps", choices = s_options, selected=sel)  
      
    })
    
    
# ----------------------------------------------------------------------------------
# ------ ROOT MAP CREATION ---------------------------------------------------------
# ----------------------------------------------------------------------------------
      
    observe({
      
      if(input$runROOTEXP == 0 | is.null(gene)){return()}
      rs.agg <- ddply(rs.melt, .(line, variable), summarize, value=mean(value))
      
  # Make the Root Map for the gene
      
      # Modifyt the SVG for each tissue
      root2 <- xmlParse("./www/root.svg")
      temp <- gene[gene$Gene_ID == input$ref_genes,]
      temp <- melt(temp, id.vars = c("Gene_ID", "prediction"))

      co <- range01(temp$value)
      temp$col <- round((co + 1)*100)
      temp$col <- (temp$col - min(temp$col))+1
      my.col <- colorRampPalette(pal)(max(temp$col))
      
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
      saveXML(root2, file='./www/root_gene_1.svg')      
      root2.bm <- rsvg("./www/root_gene_1.svg"); writePNG(root2.bm, "./www/root_gene_1.png")
      
  # Make the Root Map for the corresponding reporter
      
      # Modifyt the SVG for each tissue
      root2 <- xmlParse("./www/root.svg")
      temp <- gene[gene$Gene_ID == input$ref_genes,]
      temp <- melt(temp, id.vars = c("Gene_ID", "prediction"))
      temp <- rs.agg[rs.agg$line == temp$prediction[1] ,]
      co <- range01(temp$value)
      temp$col <- round((co + 1)*100)
      temp$col <- (temp$col - min(temp$col))+1
      my.col.1 <- colorRampPalette(pal)(max(temp$col))
      root1 <- xmlParse("./www/root.svg")
      
      co <- range01(temp$value)
      temp$col <- round((co + 1)*100)
      temp$col <- (temp$col - min(temp$col))+1
      my.col <- colorRampPalette(pal)(max(temp$col))
      
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
      saveXML(root2, file='./www/root_pred_1.svg')      
      root2.bm <- rsvg("./www/root_pred_1.svg"); writePNG(root2.bm, "./www/root_pred_1.png")      
      
      
      # Store the results
      
      root_pred_1 <<- './www/root_pred_1.png'
      root_gene_1 <<- './www/root_gene_1.png'      
      
    })
    
    
    
    # ----------------------------------------------------------------------------------
    # ------ ROOT MAP CREATION ---------------------------------------------------------
    # ----------------------------------------------------------------------------------
    
    
    observe({
        
      if(input$runROOTEXP == 0){return()}
      
      message("Creating the root maps")
      
      rs.agg <- ddply(rs.melt, .(line, variable), summarize, value=mean(value))
      rs.agg$type <- F
   
    
       
    # Make the Root Map for the wild type

      # Modifyt the SVG for each tissue
      temp <- rs.agg[rs.agg$line == input$ref_reps ,]
      co <- range01(temp$value)
      temp$col <- round((co + 1)*100)
      temp$col <- (temp$col - min(temp$col))+1
      my.col.1 <- colorRampPalette(pal)(max(temp$col))
      root1 <- xmlParse("./www/root.svg")
      i <- 1
      for(n in tissues){
        node = xpathApply(root1, paste0("//*[@id='",n,"']/*"))
        sapply(node, function(x) {
          oldstyle <- xmlAttrs(x)
          removeAttributes(x, "style")
          xmlAttrs(x)<-c(gsub("fill:#fff", paste0("fill:",my.col.1[temp$col[temp$variable == n]]),
                oldstyle['style']))    
        })
        i <- i+1
      }
      # Save the new file
      saveXML(root1, file='./www/root_1.svg')
      
      
    # Make the Root Map for the mutant
      
      # Modifyt the SVG for each tissue
      root2 <- xmlParse("./www/root.svg")
      temp <- rs.agg[rs.agg$line == input$to_plot ,]
      co <- range01(temp$value)
      temp$col <- round((co + 1)*100)
      temp$col <- (temp$col - min(temp$col))+1
      my.col.1 <- colorRampPalette(pal)(max(temp$col))
      for(n in tissues){
        node = xpathApply(root2, paste0("//*[@id='",n,"']/*"))
        sapply(node, function(x) {
          oldstyle <- xmlAttrs(x)
          removeAttributes(x, "style")
          xmlAttrs(x)<-c(gsub("fill:#fff", paste0("fill:",my.col.1[temp$col[temp$variable == n]]),
                              oldstyle['style']))    
        })
      }
      # Save the new file
      saveXML(root2, file='./www/root_2.svg')
      
      
    # Make the Root Map for the mutant-wildtype
      
      
      # Get the difference with the wild type
      rs.agg$wt <- rep(rs.agg$value[rs.agg$line == input$ref_reps], length(p.list))
      rs.agg$diff <- rs.agg$value - rs.agg$wt
      rs.agg$x <- as.numeric(rs.agg$variable)   
      
      # Create a color palette for the gene expression
      temp <- rs.agg[rs.agg$line == input$to_plot,]
      temp$pvalue <- rs.aov$pvalue[(rs.aov$genotype_1 == input$to_plot & rs.aov$genotype_2 == input$ref_reps) | 
                                  (rs.aov$genotype_1 == input$ref_reps & rs.aov$genotype_2 == input$to_plot)]
      temp$sig <- as.numeric(temp$pvalue) < 0.05
      temp <- temp[temp$sig == T,]
      tiss.diff <<- unique(as.character(temp$variable))
      
      if(nrow(temp) > 0){
        co1 <- range01(temp$diff)
        temp$col <- round((co1 + 1)*100)
        temp$col <- (temp$col - min(temp$col))+1
        my.col.2 <<- colorRampPalette(pal.div)(max(temp$col))
        
        
        # Modifyt the SVG for each tissue
        root_diff <- xmlParse("./www/root.svg")
        for(n in tiss.diff){
          node = xpathApply(root_diff, paste0("//*[@id='",n,"']/*"))
          sapply(node, function(x) {
            oldstyle <- xmlAttrs(x)
            removeAttributes(x, "style")
            xmlAttrs(x)<-c(gsub("fill:#fff", paste0("fill:",my.col.2[temp$col[temp$variable == n]]),
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
      
      
      # Store the results
      
      nameWt <<- input$ref_reps
      nameMt <<- input$to_plot
      root_1 <<- './www/root_1.png'
      root_2 <<- './www/root_2.png'
      root_diff <<- './www/root_diff.png'
      rs.agg <<- rs.agg
      
      
    })
 
    
    
    
    
# ----------------------------------------------------------------------------------
# ------ PLOT THE DATA     ---------------------------------------------------------
# ----------------------------------------------------------------------------------
    
    ## BARPLOT #############################
    
    barplot_comp_1 <- function(){
      message("-------------------------------------------")
      if(input$runROOTEXP == 0){return()}
      
      temp <- melt(gene[gene$Gene_ID == input$ref_genes,], id.vars = c("Gene_ID", "prediction"))
      temp2 <- rs.agg.short[rs.agg.short$line == temp$prediction[1],]
      temp <- temp[,-2]
      colnames(temp) <- colnames(temp2)
      temp <- rbind(temp, temp2)
      
      plot1 <- ggplot(temp, aes(variable, value, fill=line)) + 
        geom_bar(stat = "identity", position=position_dodge(width=0.9), width=0.8) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") + ylab("Relative expression value")
      #p <- ggplotly(plot1)
      plot1

    }
    output$barplot_comp_1 <- renderPlot({
      print(barplot_comp_1())
    })  
    
    ## BARPLOT #############################
    
    barplot_comp <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      temp <- rs.agg.short[rs.agg.short$line == input$ref_reps | rs.agg.short$line == input$to_plot,]
      
      plot1 <- ggplot(temp, aes(variable, value, fill=line)) + 
        geom_bar(stat = "identity", position=position_dodge(width=0.9), width=0.8) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        xlab("") + ylab("Relative expression value")
      #p <- ggplotly(plot1)
      plot1

    }
    output$barplot_comp <- renderPlot({
      print(barplot_comp())
    })      
    
  ## HEATMAP #############################
    
    heatmap <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      temp <- rs.maov
      temp[temp > 0.05] <- 0.06
      temp[temp == 0] <- 0.06
 
      dat <- as.data.frame(temp)
      dat$line_1 <- rownames(dat)
      dat <- melt(dat, id.vars = c("line_1"))
      dat$line_2 <- dat$variable
      dat$p_value <- dat$value
      ## Example data
      plot1 <- ggplot(dat, aes(line_1, line_2, z= p_value)) + geom_tile(aes(fill = p_value)) + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        scale_fill_gradient(low="blue", high="white", space="Lab")   +
        xlab("") + ylab("") + 
        coord_fixed()
      
      plot1
    }
    output$heatmap <- renderPlot({
      print(heatmap())
    })   
 
    
    
  ## ROOT 1 #############################
    
    output$nameWt <- renderText({ 
      if(input$runROOTEXP == 0){return()}
      nameWt  
    })  
    
    output$root_1 <- renderImage({
      if(input$runROOTEXP == 0){
        filename <- normalizePath(file.path("./www/root.png"))
        return(list(src = filename, width=150))
      }      
      filename <- normalizePath(file.path(root_1))
      list(src = filename, width=150)     
    }, deleteFile = FALSE
    )
    
    
    
  ## ROOT 2 #############################
    
    output$nameMt <- renderText({
      if(input$runROOTEXP == 0){return()}
      nameMt  
    })  
    
    output$root_2 <- renderImage({
      if(input$runROOTEXP == 0){
        filename <- normalizePath(file.path("./www/root.png"))
        return(list(src = filename, width=150))
      }
      filename <- normalizePath(file.path(root_2))
      list(src = filename, width=150)     
    }, deleteFile = FALSE
    )
    
    
  ## ROOT DIFFERENCE #############################
    
    output$nameDiff <- renderText({  
      if(input$runROOTEXP == 0){return()}
      paste0(nameMt, " - ", nameWt)  
    }) 
    
    output$root_diff <- renderImage({
      if(input$runROOTEXP == 0){
        filename <- normalizePath(file.path("./www/root.png"))
        return(list(src = filename, width=150))
      }
      filename <- normalizePath(file.path(root_diff))
      list(src = filename, width=150)     
    }, deleteFile = FALSE
    )  
    
    ## ROOT GENE #############################
    
    output$name_gene_1 <- renderText({
      if(input$runROOTEXP == 0){return()}
      input$ref_genes
    })  
    
    output$root_gene_1 <- renderImage({
      if(input$runROOTEXP == 0){
        filename <- normalizePath(file.path("./www/root.png"))
        return(list(src = filename, width=150))
      }
      filename <- normalizePath(file.path(root_gene_1))
      list(src = filename, width=150)     
    }, deleteFile = FALSE
    )    
    
    ## ROOT GENE PREDICT #############################
    
    output$name_pred_1 <- renderText({
      if(input$runROOTEXP == 0){return()}
      as.character(gene$prediction[gene$Gene_ID == input$ref_genes])
    })  
    
    output$root_pred_1 <- renderImage({
      if(input$runROOTEXP == 0){
        filename <- normalizePath(file.path("./www/root.png"))
        return(list(src = filename, width=150))
      }
      filename <- normalizePath(file.path(root_pred_1))
      list(src = filename, width=150)     
    }, deleteFile = FALSE
    )    
    
 
  ## FLUORESCENCE DISTRIBUTIONS #############################
  
    histoPlot <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      temp <- melt(rs, .(line, root))
      message(str(temp))
      plot1 <- ggplot(temp, aes(x=value, fill=line, colour=line)) +
        geom_density(alpha=0.4, lwd=0.8, adjust=0.5) +
        theme_bw() 
      plot1
    }
    
    output$histoPlot <- renderPlot({
      print(histoPlot())
    })
   
     
    output$downloadPlot <- downloadHandler(
      filename = "distribution.png",
      content = function(file) {
        png(file, width = 400, height=300)
        histoPlot()
        dev.off()
      }
    )
  
    ## FLUORESCENCE DISTRIBUTIONS #############################
    
    ldaPlot <- function(){
      
      if(input$runROOTEXP == 0){return()}
      
      temp <- rs
      temp$col <- 0
      temp$col[temp$line == input$to_plot] <- 1
      temp$col[temp$line == input$ref_reps] <- 2
      temp$line[temp$line != input$to_plot & temp$line != input$ref_reps] <- " "
      
      plot2 <- ggplot(temp[temp$col>0,], aes(LD1, LD2, colour=line)) + 
        geom_point(data=temp[temp$col==0,], aes(LD1, LD2), colour="grey", size=2) +
        geom_point(size=4) +
        theme_bw() + 
        scale_colour_manual(values=c("#F96F70", "#00BC47"), 
                          name="Lines",
                          labels=c(input$to_plot, input$ref_reps)) +
        stat_ellipse(level = 0.6)
      plot2
    }
    
    output$ldaPlot <- renderPlot({
      if(input$runROOTEXP == 0){return()}
      print(ldaPlot())
    })
    
    
    
    ## COLOR SCALE ################################################
    
    output$scalePlot <- renderPlot({

      y <- seq(from=0, to=1, length.out = 100)
      y2 <- seq(from=0, to=1, length.out = 5)
      x <- rep(1, length(100))
      x2 <- rep(0.9, 5)
      temp.col <- data.frame(x, y)
      temp.col.2 <- data.frame(x2, y2)
      plot1 <- ggplot(temp.col, aes(x, y, colour=factor(y))) + 
        theme_void() + 
        theme(legend.position="none", text = element_text(size=30), plot.title=element_text(size=14)) + 
        geom_point(pch=15, size=4) + 
        scale_colour_manual(values=colorRampPalette(pal)(100)) +
        geom_text(data=temp.col.2, aes(x=x2, y=y2, label=y2), colour="black", size=5) + 
        xlim(c(0.8, 1.05)) + 
        ggtitle("\n\n \n \n")
      
      y <- seq(from=-1, to=1, length.out = 100)
      y2 <- seq(from=-1, to=1, length.out = 5)
      temp.col <- data.frame(x, y)
      temp.col.2 <- data.frame(x2, y2)
      plot2 <- ggplot(temp.col, aes(x, y, colour=factor(y))) + 
        theme_void() + 
        theme(legend.position="none", text = element_text(size=30), plot.title = element_text(size=15)) + 
        geom_point(pch=15, size=4) + 
        scale_colour_manual(values=colorRampPalette(pal.div)(100)) +
        geom_text(data=temp.col.2, aes(x=x2, y=y2, label=y2), colour="black", size=5) + 
        xlim(c(0.8, 1.05)) + 
        ylim(c(-1, 1)) + 
        ggtitle("\n\n \n \n")
      
      remove(temp.col, temp.col.2, x, y, x2, y2)
      pl <- grid.arrange(plot1, plot2, nrow=2)
      print(pl)
    })
    
    output$scalePlot2 <- renderPlot({
      
      y <- seq(from=0, to=1, length.out = 100)
      y2 <- seq(from=0, to=1, length.out = 5)
      x <- rep(1, length(100))
      x2 <- rep(0.9, 5)
      temp.col <- data.frame(x, y)
      temp.col.2 <- data.frame(x2, y2)
      plot1 <- ggplot(temp.col, aes(x, y, colour=factor(y))) + 
        theme_void() + 
        theme(legend.position="none", text = element_text(size=30), plot.title=element_text(size=14)) + 
        geom_point(pch=15, size=4) + 
        scale_colour_manual(values=colorRampPalette(pal)(100)) +
        geom_text(data=temp.col.2, aes(x=x2, y=y2, label=y2), colour="black", size=5) + 
        xlim(c(0.8, 1.05)) + 
        ggtitle("\n\n \n \n")
    
      remove(temp.col, temp.col.2, x, y, x2, y2)
      print(plot1)
    })
    
    
    
    
    output$line_comp_text <- renderText({ 
      if(input$runROOTEXP == 0){return()}
      sig <- rs.maov[input$ref_reps, input$to_plot]
      text <- ""
      if(sig <= 0.05){ text <- "The overall difference between lines IS statistically significant."
      }else{ text <- "The overlall difference between lines IS NOT statistically significant."}
      return(HTML(text))
    })
    
    
    output$line_comp_pval <- renderText({ 
      if(input$runROOTEXP == 0){return()}
      
      sig <- rs.maov[input$ref_reps, input$to_plot]
      text <- ""
      if(sig <= 0.05){ text <- paste0("p-value = ",round(sig, 3))
      }else{text <- paste0("p-value = ",round(sig, 3))}
      return(HTML(text))
    })    

    
    output$cell_comp <- renderText({ 
      if(input$runROOTEXP == 0){return()}
      temp <- rs.agg[rs.agg$line == input$to_plot,]
      temp$pvalue <- rs.aov$pvalue[(rs.aov$genotype_1 == input$to_plot & rs.aov$genotype_2 == input$ref_reps) | 
                                     (rs.aov$genotype_1 == input$ref_reps & rs.aov$genotype_2 == input$to_plot)]
      temp$sig <- as.numeric(temp$pvalue) < 0.05
      temp <- temp[temp$sig == T,]
      unique(as.character(temp$variable))
    })
    
    
    
    
    #------------------------------------------------------
    #------------------------------------------------------
    # TABLES
    #------------------------------------------------------
    #------------------------------------------------------ 
    
    # MAOV RESULTS 
    
    output$maov_results <- renderTable({
      if (input$runROOTEXP == 0) { return()}
      dat <- as.data.frame(rs.maov)
      dat$line_1 <- rownames(dat)
      dat <- melt(dat, id.vars = c("line_1"))
      colnames(dat) <- c("line_1", "line_2", "p_value")
      dat
    })    
    output$download_moav <- downloadHandler(
      filename = function() {"moav_results.csv"},
      content = function(file) {
        dat <- as.data.frame(rs.maov)
        dat$line_1 <- rownames(dat)
        dat <- melt(dat, id.vars = c("line_1"))
        colnames(dat) <- c("line_1", "line_2", "p_value")
        write.csv(dat, file)
      }
    )
    
    # AOV RESULTS 
    
    output$aov_results <- renderTable({
      if (input$runROOTEXP == 0) { return()}
      rs.aov
    })    
    output$download_oav <- downloadHandler(
      filename = function() {"oav_results.csv"},
      content = function(file) {
        write.csv(rs.aov, file)
      }
    )
    
    # COMPARISON RESULTS 
    
    output$comp_results <- renderTable({
      if (input$runROOTEXP == 0) { return()}
      gene
    })    
    output$download_comp <- downloadHandler(
      filename = function() {"oav_results.csv"},
      content = function(file) {
        write.csv(gene, file)
      }
    )    
    
    
    #     
    #     
    #     
    #     
    #     # Plot the different distributions
    #     diffPlot <- function(){
    #       
    #       if(input$runROOTEXP == 0){return()}
    #       
    #       rs.agg <- Results()$agg
    #       rs.agg$mean.diff <- Results()$mean.diff
    #       tissues <- Results()$tissues
    # 
    #       plot1 <- ggplot(data = rs.agg, aes(x, diff, colour=plant)) + 
    #         geom_rect(aes(xmin=min(x), xmax=max(x), 
    #                       ymin=-mean.diff, ymax=mean.diff), 
    #                   fill="lightgrey", color="lightgrey", alpha=0.5) +
    #         geom_point(data=rs.agg[rs.agg$type == 1,], 
    #                    aes(x, diff, colour=plant), size=6, pch=1) +
    #         geom_line(size=1) + 
    #         geom_point(size=3) +
    #         scale_x_discrete(breaks= c(1:length(tissues)),labels=tissues) +
    #         theme_bw()
    #       plot1
    # 
    #     }
    #     
    #     output$diffPlot <- renderPlot({
    #       if(input$runROOTEXP == 0){return()}
    #       
    #       print(diffPlot())
    #     })
    #     
    #     output$downloadPlot <- downloadHandler(
    #       filename = "difference.png",
    #       content = function(file) {
    #         png(file, width = input$plot_width, height=input$plot_height)
    #         diffPlot()
    #         dev.off()
    #       })     
    #     
    #     
    #     
    #     
    #     # Plot the different distributions
    #     diffPlotSingle <- function(){
    #       
    #       if(input$runROOTEXP == 0){return()}
    #       
    #       rs.agg <- Results()$agg
    #       rs.agg$mean.diff <- Results()$mean.diff
    #       tissues <- Results()$tissues
    #       
    #       plot1 <- ggplot(rs.agg, aes(x, diff, colour=plant)) + 
    #         geom_rect(aes(xmin=min(rs.agg$x), xmax=max(rs.agg$x), 
    #                       ymin=-mean.diff, ymax=mean.diff), 
    #                   fill="lightgrey", color="lightgrey", alpha=0.5) +
    #         geom_line(size=0.5) + 
    #         geom_line(data = rs.agg[rs.agg$plant == input$to_plot,], aes(x, diff), colour="#FF7C00", size=1) + 
    #         geom_point(data=rs.agg[rs.agg$type == 1 & rs.agg$plant == input$to_plot,], aes(x, diff), colour="#FF7C00",size=6, pch=1) +
    #         geom_point(size=2) +
    #         geom_point(data=rs.agg[rs.agg$plant == input$to_plot,], aes(x, diff), colour="#FF7C00", size=3) +
    #         scale_x_discrete(breaks= c(1:length(tissues)),labels=tissues) +
    #         scale_colour_grey()+
    #         theme_bw() +
    #         theme(axis.text.x = element_text(angle = 90, hjust = 1))
    #       
    #       plot1
    #       
    #     }
    #     output$diffPlotSingle <- renderPlot({
    #       if(input$runROOTEXP == 0){return()}
    #       
    #       print(diffPlotSingle())
    #     })    
    #     
    #       
    
    #       temp2 <- Results()$norm
    #       plot2 <- ggplot(temp2, aes(LD1, LD2, colour=line)) + 
    #         geom_point(size=4) +
    #         theme_bw() + 
    #         stat_ellipse()
    

})
