
#
# Copyright © 2017, Université catholique de Louvain
# All rights reserved.
# 
# Copyright © 2017 Forschungszentrum Jülich GmbH
# All rights reserved.
# 
# Developers: Guillaume Lobet
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted under the GNU General Public License v3 and provided that the following conditions are met:
#   
# 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
# 
# Disclaimer
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# You should have received the GNU GENERAL PUBLIC LICENSE v3 with this file in license.txt but can also be found at http://www.gnu.org/licenses/gpl-3.0.en.html
# 
# NOTE: The GPL.v3 license requires that all derivative work is distributed under the same license. That means that if you use this source code in any other program, you can only distribute that program with the full source code included and licensed under a GPL license.


# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com

library(shiny)

# Source the external file contaiing the plot functions
source("www/plot_functions.R")
options(shiny.maxRequestSize=30*1024^2) 
gene <<- NULL



shinyServer(
  function(input, output, clientData, session) {

    load("www/example/root.RData")
    root$value[root$tissue != "borders"] = 1
    root$value[root$tissue == "borders"] = NA
    rs <- reactiveValues(reporter = NULL, 
                         gene = NULL,
                         lda.fit = NULL,
                         rep.melt = NULL,
                         tissues = NULL,
                         rep.aov = NULL,
                         rep.aov.gene = NULL,
                         rep.maov = NULL,
                         p.list= NULL,
                         gene.list = NULL,
                         rep.agg.short = NULL,
                         gene.dist = NULL,
                         root = root,
                        dts = NULL)
    
    
    observe({
      
      fileName <- 'www/datasets.txt'
      dts <- read.table(fileName, sep="\t", stringsAsFactors = F, header = T)
      rs$dts <- dts

      fileName <- 'www/datasets/reporters/'
      flist <- unique(gsub(".csv", "", list.files(fileName)))
      fl <- dts$name[dts$id %in% flist]
      updateSelectInput(session, "reporters", choices = fl)  
      
      fileName <- 'www/datasets/microarrays/'
      flist <- unique(gsub(".csv", "", list.files(fileName)))
      fl <- dts$name[dts$id %in% flist]
      updateSelectInput(session, "microarrays", choices = fl)  
    })         
    
    
    #------------------------------------------------------
    #------------------------------------------------------
    #         LOAD THE USER DATA
    #------------------------------------------------------
    #------------------------------------------------------
     
    observeEvent(input$load_data, {
            

      withProgress(message = 'WORKING:', value = 0, {
      
        
        incProgress(1/4, detail = "Loading the data")
        

        # if(input$use_example){ # Use the exmple dataset
         # gene <- fread("www/datasets/microarrays/all_genes.csv")
         # temp <- read_rsml("www/datasets/reporters/experession_full.rsml")
         # write.csv(temp, "www/datasets/reporters/experession_full.csv", row.names = F)
         # temp2 <- read.csv("www/datasets/reporters/experession_full.csv")
        #   # load("www/root.RData")
        #   # load("www/example.RData")
        #   for(l in list.files("www/example/")) load(paste0("www/example/", l))
        #   #rs <- rs
        # }else{
          
          # Load two datafiles, for the gene and the reporters

          if(input$use_example){
            dat <- rs$dts$id[rs$dts$name == input$microarrays]
            gene <- fread(paste0('www/datasets/microarrays/', dat,".csv"))
            dat <- rs$dts$id[rs$dts$name == input$reporters]
            temp <- read.csv(paste0('www/datasets/reporters/', dat,".csv"))
          }else{
            inGene <- input$gene_file
            inReporter <- input$rep_file

            if (is.null(inReporter)) return(NULL)
            if(!is.null(inGene)){
              gene <- fread(inGene$datapath)   
            }
            temp <- read_rsml(inReporter$datapath)   
          }
          
          
          # Average the data by line, root, cell type
          if(input$method == "Mean") mean_data <- ddply(temp, .(line, root, cell_type), summarise, value=mean(value))
          if(input$method == "Median") mean_data <- ddply(temp, .(line, root, cell_type), summarise, value=median(value))
          if(input$method == "Min") mean_data <- ddply(temp, .(line, root, cell_type), summarise, value=min(value))
          if(input$method == "Max") mean_data <- ddply(temp, .(line, root, cell_type), summarise, value=max(value))
          mean_data <- mean_data[!is.na(mean_data$value),]
          
          # Average the data by line, cell type
          if(input$method == "Mean") mean_data_2 <- ddply(temp, .(line, cell_type), summarise, value=mean(value))
          if(input$method == "Median") mean_data_2 <- ddply(temp, .(line, cell_type), summarise, value=median(value))
          if(input$method == "Min") mean_data_2 <- ddply(temp, .(line, cell_type), summarise, value=min(value))
          if(input$method == "Max") mean_data_2 <- ddply(temp, .(line, cell_type), summarise, value=max(value))
          mean_data_2 <- mean_data_2[!is.na(mean_data$value),]          
          
          # Reshape the data to have them in the proper form for the analysis
          reporter <- dcast(mean_data, line + root ~ cell_type)
          for(cn in colnames(reporter)){
            reporter[[cn]][is.na(reporter[[cn]])] <- 0
          }
          
          rep.mean <- dcast(mean_data_2, line ~ cell_type) # For the display
        
        
        
        #------------------------------------------------------
        # PROCESS THE DATA
        #------------------------------------------------------
  
        p.list <- as.character(unique(reporter$line))
        gene.list <- as.character(unique(gene$Gene_ID))
        tissues <- colnames(reporter)[-c(1,2)]
        
        # Aggregate the data   
        rep.melt <- melt(reporter, id.vars =c("line", "root"))
        # gene.melt <- melt(gene, id.vars =c("Gene_ID"))
        rep.agg.short <- ddply(rep.melt, .(line, variable), summarize, value=mean(value))
        
        # ----------------------------------------------------------------------------------
        # ----- PAIRWISE ANOVA COMPARISONS AND MANOVA ANALYSIS   ---------------------------
        # ----------------------------------------------------------------------------------            
        
        incProgress(1/4, detail = "Analysing the data")
        
          # Create a table that will contain the anova results
          l1 <- length(p.list)
          l2 <- l1^2
          l3 <- (((l1 * l1) - l1) / 2) * length(tissues)
          
          aov.results <<- data.frame(genotype_1 = character(l3), 
                                     genotype_2=character(l3), 
                                     tissue = character(l3), 
                                     pvalue = numeric(l3), 
                                     stringsAsFactors = F)
          
          maov.results <<- matrix(NA, ncol = l1, nrow = l1) 
          colnames(maov.results) <- p.list 
          rownames(maov.results) <- p.list
          
          i <- 1
          k <- 1
          for(p in p.list){
            j <- 1
            for(p1 in p.list){
              if(p != p1 & is.na(maov.results[j,i])){
                
                # MANOVA analysis to compare the plants
                temp <- reporter[reporter$line == p | reporter$line == p1,]
                
                # Get the ines sleected by the user
                ts <- c("columella","cortex","endodermis","epidermis","lateralrootcap","QC","stele")
                # ts <- input$type_to_analyse
                # This is a very ugly way to do this; 
                for(t in c(1:length(ts))) temp[[paste0("V",t)]] <- temp[[ts[t]]]
                if(length(ts) == 1) fit <- manova(cbind(V1) ~ line,  temp)
                if(length(ts) == 2) fit <- manova(cbind(V1,V2) ~ line,  temp)
                if(length(ts) == 3) fit <- manova(cbind(V1,V2,V3) ~ line,  temp)
                if(length(ts) == 4) fit <- manova(cbind(V1,V2,V3,V4) ~ line,  temp)
                if(length(ts) == 5) fit <- manova(cbind(V1,V2,V3,V4,V5) ~ line,  temp)
                if(length(ts) == 6) fit <- manova(cbind(V1,V2,V3,V4,V5,V6) ~ line,  temp)
                if(length(ts) == 7) fit <- manova(cbind(V1,V2,V3,V4,V5,V6,V7) ~ line,  temp)
                # fit <- manova(cbind(lateralrootcap, columella, QC, epidermis, cortex, endodermis, stele) ~ line,  temp) # Stele was removed

                maov.results[j,i] <- tryCatch({
                  round(summary(fit, test = "Wilks", tol=0)$stats[1,6], 5)
                },warning = function(w) {
                }, error = function(e) {
                  -1
                })
                maov.results[i,j] <- maov.results[j,i]

                # Reponse for each tissue
                for(ti in tissues){
                  temp <- rep.melt[rep.melt$line == p | rep.melt$line == p1,]
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
        
          
          # ----------------------------------------------------------------------------------
          # ----- PAIRWISE DISTANCE COMPARISONS BETWEEN GENES AND REPORTERS   ---------------------------
          # ----------------------------------------------------------------------------------            
          
          incProgress(1/5, detail = "Matching the genes")
          
          # # Create a table that will contain the anova results
          # l1 <- length(p.list)
          # l2 <- length(gene.list)
          # 
          # dist.results <<- matrix(0, ncol = l1, nrow = l2) 
          # colnames(dist.results) <- p.list 
          # rownames(dist.results) <- gene.list
          # 
          # i <- 1
          # for(p in p.list){
          #   j <- 1
          #   for(p1 in gene.list){
          #     if(dist.results[j,i] == 0){
          #       
          #       # MANOVA analysis to compare the plants
          #       temp <- rep.agg.short[rep.agg.short$line == p,]
          #       temp2 <- gene.melt[gene.melt$Gene_ID == p1,]
          #       
          #       # Get the ines sleected by the user
          #       #ts <- c("columella","cortex","endodermis","epidermis")#,"lateralrootcap","QC","stele")#input$type_to_analys
          #       ts <- input$type_to_analyse
          #       
          #       temp <- temp[temp$variable %in% ts,]
          #       temp2 <- temp2[temp2$variable %in% ts,]
          #       
          #       dist.results[j,i] <- dist(rbind(temp$value, temp2$value))
          #       #dist.results[i,j] <- dist.results[j,i]
          #       
          #     }
          #     j <- j+1
          #   }
          #   i <- i+1
          # }
          # 
          # dist <- data.frame(dist.results)
          # dist$gene <- rownames(dist)
          # dist <- melt(dist, id.vars=c("gene"))
          # 
          # for(g in gene.list){
          #   temp <- dist[dist$gene == g,]
          #   temp <- temp[temp$value  == min(temp$value),]
          #   gene$match[gene$Gene_ID == g] <- as.character(temp$variable[1])
          #   gene$distance[gene$Gene_ID == g] <- temp$value[1]
          # }
          # 
          # dist <- data.frame(t(dist.results))
          # dist$rep <- rownames(dist)
          # dist <- melt(dist, id.vars=c("rep"))
          aov.results.gene <- NULL
          for(p in p.list){
            ge <- strsplit(p, "_")[[1]]
            ge <- ge[grepl("AT[1-9]G", ge)] # find the corresping gene
            ge <- gsub("AT", "At", ge)
            ge <- gsub("G", "g", ge)
            # temp <- temp[temp$value  == min(temp$value),]
            if(length(ge) > 0){
              rep.mean$match[rep.mean$line == p] <- ge #as.character(temp$variable[1])
              # rep.mean$distance[rep.mean$line == p] <- temp$value[1]
              rep.agg.short$match[rep.agg.short$line == p] <- ge #as.character(temp$variable[1])
              rep.melt$match[rep.melt$line == p] <- ge #as.character(temp$variable[1])
              # rep.agg.short$distance[rep.agg.short$line == p] <- temp$value[1]
              
              temp <- rep.melt[rep.melt$line == p,]
              temp2 <- gene[gene$Gene_ID == ge,]
              colnames(temp2) <- c("line", "variable", "value")
              for(t in tissues){
                tmp <- rbind(temp[temp$variable == t,c("line", "value")],temp2[temp2$variable == t,c("line", "value")])
                fit <- aov(value ~ line, data=tmp)
                aov.results.gene <- rbind(aov.results.gene, data.frame(line=p, gene=ge, tissue=t, pvalue = round(summary(fit)[[1]][1,5], 5)))
              }
            }
          }
          
          
          
        # -------------------------------------------------------
        # ----- LDA ANALYSIS ------------------------------------
        # -------------------------------------------------------            
        
        incProgress(1/5, detail = "Performing the PCA")
        
          # Make the LDA analysis on the reporter dataset.
          # ts <- c("line", input$type_to_analyse)
          temp <- reporter[complete.cases(reporter),]
          lines <- temp$line
          temp <- temp[,input$type_to_analyse]
          
          
          
          pca <- prcomp(temp, retx = T, scale=T)  # Make the PCA
          pca.results <- cbind(line=lines, data.frame(pca$x)[,])
          
          # reporter <- 
          # 
          # 
          # temp$line <- factor(temp$line) # To avoid empty group if they occur
          # fit <- lda(line ~ ., data=temp, CV=F)
          # 
          # 
          # 
          # 
          # 
          # # Get the accuracy of the prediction
          # fit.p <- predict(fit, newdata = temp[,-1])
          # ct <- table(temp$line, fit.p$class)
          # #diag(prop.table(ct, 1))
          # #sum(diag(prop.table(ct)))
          # 
          # #gene <- read.csv("~/Desktop/GeneExpression.txt", sep="\t")
          # if(!is.null(gene)){
          #   prediction <- predict(fit, newdata=gene[,-1])$class
          #   gene <- cbind(gene, prediction)
          # }
          # 
          # reporter <- cbind(reporter, fit.p$x)
        # }
        
          
          
        # save(gene, file="www/example/gene.RData")
        # save(fit, file="www/example/fit.RData")
        # save(reporter, file="www/example/rep.RData")
        # save(rep.melt, file="www/example/repmelt.RData")
        # save(tissues, file="www/example/tissues.RData")
        # save(aov.results, file="www/example/aov.RData")
        # save(maov.results, file="www/example/maov.RData")
        # save(p.list, file="www/example/plist.RData")
        # save(rep.agg.short, file="www/example/repaggshort.RData")
        
        
          rs$reporter <- reporter
          rs$rep.mean <- rep.mean
          rs$gene <- gene
          #rs$lda.fit <- fit
          rs$lda.fit <- pca.results
          rs$reporter <- reporter
          rs$rep.melt <- rep.melt
          rs$tissues <- tissues
          rs$rep.aov <- aov.results
          rs$rep.aov.gene <- aov.results.gene
          rs$rep.maov <- maov.results
          rs$p.list <- p.list
          rs$rep.agg.short <- rep.agg.short
          # rs$gene.dist <- dist.results

      #------------------------------------------------------
      # ----- UPDATE THE UI DATA ----------------------------
      #------------------------------------------------------

      incProgress(1/4, detail = "Updating the UI")

        # reporter <- rs$reporter        
        reps <- na.omit(unique(factor(reporter$line)))
        
        # Genotype list
        s_options <- list()
        for(r in reps) s_options[[r]] <- r
        updateSelectInput(session, "ref_reps", choices = s_options)  
        updateSelectInput(session, "ref_reps_2", choices = s_options)  
        
        # Genotype check box
        updateSelectInput(session, "to_plot", choices = s_options, selected=s_options[2])  
        
        # gene <- rs$gene
        # # Genes
        # if(!is.null(gene)) {
        #   gens <- na.omit(factor(unique(gene$Gene_ID)))
        #   g_options <- list()
        #   for(g in gens) g_options[[g]] <- g
        #   updateSelectInput(session, "ref_genes", choices = g_options)        
        # }
      })
    })
    
    
    
  # ----------------------------------------------------------------------------------
  # ------ UPDATE UI TO PREVENT IDENTIQUE SELECTIONS           ------------------
  # ----------------------------------------------------------------------------------
    
    observe({
      
      if(is.null(rs$reporter) || grepl("Load", input$to_plot)){return()}
      
      sel <- input$to_plot
      reps <- na.omit(unique(factor(rs$reporter$line)))
      reps1 <- reps[-match(input$ref_reps, reps)]

      # Genotype check box
      s_options <- list()
      for(r in reps1) s_options[[r]] <- r
      updateSelectInput(session, "to_plot", choices = s_options, selected=sel)  
    })
    
    observe({
      
      ct_options <- list()
      sel <- input$type_to_analyse
      if(length(sel) == 0) sel = cell_types
      for(ct in cell_types) ct_options[[ct]] <- ct
      updateSelectInput(session, "type_to_analyse", choices = ct_options, selected=sel) 
    })
    
    observe({
      
      if(is.null(rs$reporter) || grepl("Load", input$ref_reps)){return()}
      
      sel <- input$ref_reps
      reps <- na.omit(unique(factor(rs$reporter$line)))
      reps2 <- reps[-match(input$to_plot, reps)]
      
      # Genotype list
      s_options <- list()
      for(r in reps2) s_options[[r]] <- r
      updateSelectInput(session, "ref_reps", choices = s_options, selected=sel)  
      
    })
    
    
    observe({
      
      if(is.null(rs$reporter) || grepl("Load", input$ref_reps_2)){return()}
      
      sel <- input$ref_reps_2
      reps <- na.omit(unique(factor(rs$reporter$line)))
      
      # Genotype list
      s_options <- list()
      for(r in reps) s_options[[r]] <- r
      updateSelectInput(session, "ref_reps_2", choices = s_options, selected=sel)  
      
    })    
    
    
# ----------------------------------------------------------------------------------
# ------ PLOT THE DATA     ---------------------------------------------------------
# ----------------------------------------------------------------------------------
    

## BARPLOT #############################
    
    output$barplot_comp_1 <- renderPlot({
      if(is.null(rs$gene)){return()}
      print(barplot_comp_1(reps = input$ref_reps_2,
                           rep.agg = rs$rep.melt, 
                           gene = rs$gene
                            ))
    })  
 
    
## BARPLOT #############################
    
    output$barplot_comp <- renderPlot({
      if(is.null(rs$reporter)){return()}
      print(barplot_comp(input$ref_reps, input$to_plot, rs$rep.melt))
    })   

## PLOT THE ROOT #############################
    
    output$plotRootKey <- renderPlot({
      print(plotRootKey(rs$root[rs$root$id == 1,]))
    }) 

    
## PLOT THE ROOT WITH THE REPORTER DATA  #############################
    
    output$plotRoot <- renderPlot({
      if(is.null(rs$reporter)){return()}
      print(plotRootReporters(reps = input$ref_reps, 
                              to_plot = input$to_plot, 
                              rep.aov = rs$rep.aov, 
                              root = rs$root, 
                              rep.melt = rs$rep.melt, 
                              rep.agg.short = rs$rep.agg.short,
                              sig = rs$rep.maov[input$ref_reps, input$to_plot] < 0.05,
                              input$show_diff))
    }) 
    
## PLOT THE ROOT WITH THE GENE DATA  #############################
    
    output$plotRootGene <- renderPlot({
      if(is.null(rs$gene)){return()}
      print(plotRootGene(reps = input$ref_reps_2, 
                         root = rs$root, 
                         gene = rs$gene,
                         rep.agg.short = rs$rep.agg.short
                         ))
    })     
 
    
## PLOT THE LDA PLOT FOR THE REPORTERS  #############################
    
    output$ldaPlot <- renderPlotly({
      if(is.null(rs$reporter) || grepl("Load", input$to_plot)){return()}
      print(ldaPlot(reps = input$ref_reps, to_plot = input$to_plot, rep = rs$lda.fit))
    })
    
    
## PLOT THE HEATMAP OF THE MANOVA ANALYSIS  
    
    output$heatmap <- renderPlot({
      if(is.null(rs$rep.maov)){return()}
      print(heatmap(rs$rep.maov))
    })   
 
    
## PLOT THE HEATMAP OF THE DISTANCE ANALYSIS  
    
    output$heatmap_dist <- renderPlot({
      if(is.null(rs$gene.dist)){return()}
      print(heatmap_dist(rs$gene.dist))
    })   
    
    
    
# ----------------------------------------------------------------------------------
# ------ UPDATE THE TEXT FIELDS     ---------------------------------------------------------
# ----------------------------------------------------------------------------------
    
    
    output$name_gene_1 <- renderText({
      if(is.null(rs$reporter)){return()}
      input$ref_genes
    })  
    
    
    output$name_pred_1 <- renderText({
      if(is.null(rs$reporter) || grepl("Load", input$ref_genes)){return()}
      gene <- rs$gene
      as.character(gene$prediction[gene$Gene_ID == input$ref_genes])
    })  
  
    
    output$line_comp_text <- renderText({ 
      if(is.null(rs$rep.maov) || grepl("Load", input$to_plot)){return()}
      sig <- rs$rep.maov[input$ref_reps, input$to_plot]
      text <- ""
      if(!is.na(sig)){
        if(sig <= 0.05){ text <- "The overall difference between lines IS statistically significant."
        }else{ text <- "The overlall difference between lines IS NOT statistically significant."}
      }
      return(HTML(text))
    })
    
    
    output$line_comp_pval <- renderText({ 
      if(is.null(rs$rep.maov) || grepl("Load", input$to_plot)){return()}
      sig <- rs$rep.maov[input$ref_reps, input$to_plot]
      text <- ""
      if(!is.na(sig)){
        if(sig <= 0.05){ text <- paste0("p-value = ",round(sig, 3))
        }else{text <- paste0("p-value = ",round(sig, 3))}
      }
      return(HTML(text))
    })    

    
    output$cell_comp <- renderText({ 
      if(is.null(rs$rep.agg.short)){return()}
      
      temp <- rs$rep.agg.short[rs$rep.agg.short$line == input$to_plot,]
      
      rep.aov <- rs$rep.aov
      pvals <- rep.aov[(rep.aov$genotype_1 == input$to_plot & rep.aov$genotype_2 == input$ref_reps) | 
                                (rep.aov$genotype_1 == input$ref_reps & rep.aov$genotype_2 == input$to_plot),c("tissue", "pvalue")]
      temp <- merge(temp, pvals, by.x="variable", by.y = "tissue")
      temp$sig <- as.numeric(temp$pvalue) < 0.05
      temp <- temp[temp$sig == T,]
      
      paste(temp$variable, "( p-val:", round(as.numeric(temp$pvalue), 4), ")", collapse=" || ")
    })
    
    
    
    output$gene_comp <- renderText({ 
      if(is.null(rs$rep.melt)){return()}
      
      temp <- rs$rep.agg.short[rs$rep.agg.short$line == input$ref_reps_2,]
      
      rep.aov.gene <- rs$rep.aov.gene
      pvals <- rep.aov.gene[rep.aov.gene$line == input$ref_reps_2,c("tissue", "pvalue")]
      temp <- merge(temp, pvals, by.x="variable", by.y = "tissue")
      temp$sig <- as.numeric(temp$pvalue) < 0.05
      temp <- temp[temp$sig == T,]
      
      paste(temp$variable, "( p-val:", round(as.numeric(temp$pvalue), 4), ")", collapse=" || ")
    })    
    
    
    
    
    output$littTitle <- renderUI( {
      if(is.null(rs$dts)){return()}
      strong(rs$dts$title[rs$dts$name == input$reporters])
    }) 
    
    output$littAuth <- renderUI( {
      if(is.null(rs$dts)){return()}
      rs$dts$author[rs$dts$name == input$reporters]
    }) 
    
    output$littRef <- renderUI( {
      if(is.null(rs$dts)){return()}
      litt <- rs$dts[rs$dts$name == input$reporters,]
      paste0(litt$journal, ", ", litt$volume, ", ", litt$pages, ", ", litt$year)
    })     
    
    output$doi <- renderUI( {
      if(is.null(rs$dts)){return()}
      litt <- rs$dts
      link <- paste0("http://dx.doi.org/", litt$doi[litt$name == input$reporters])
      a("View paper", href=link, target="_blank")
    }) 
    
    
    
    output$littTitle1 <- renderUI( {
      if(is.null(rs$dts)){return()}
      strong(rs$dts$title[rs$dts$name == input$microarrays])
    }) 
    
    output$littAuth1 <- renderUI( {
      if(is.null(rs$dts)){return()}
      rs$dts$author[rs$dts$name == input$microarrays]
    }) 
    
    output$littRef1 <- renderUI( {
      if(is.null(rs$dts)){return()}
      litt <- rs$dts[rs$dts$name == input$microarrays,]
      paste0(litt$journal, ", ", litt$volume, ", ", litt$pages, ", ", litt$year)
    })     
    
    output$doi1 <- renderUI( {
      if(is.null(rs$dts)){return()}
      litt <- rs$dts
      link <- paste0("http://dx.doi.org/", litt$doi[litt$name == input$microarrays])
      a("View paper", href=link, target="_blank")
    })     
    
#------------------------------------------------------
#------------------------------------------------------
# TABLES
#------------------------------------------------------
#------------------------------------------------------ 
    
    # MAOV RESULTS 
    
    output$maov_results <- renderTable({
      if(is.null(rs$rep.maov)){return()}
      dat <- as.data.frame(rs$rep.maov)
      dat$line_1 <- rownames(dat)
      dat <- melt(dat, id.vars = c("line_1"))
      colnames(dat) <- c("line_1", "line_2", "p_value")
      dat[dat$line_1 == input$to_plot,]
    })   
    
    output$download_moav <- downloadHandler(
      filename = function() {"moav_results.csv"},
      content = function(file) {
        dat <- as.data.frame(rs$rep.maov)
        dat$line_1 <- rownames(dat)
        dat <- melt(dat, id.vars = c("line_1"))
        colnames(dat) <- c("line_1", "line_2", "p_value")
        write.csv(dat[dat$line_1 == input$to_plot,], file)
      }
    )
    
    output$maov_results_all <- renderTable({
      if(is.null(rs$rep.maov)){return()}
      dat <- as.data.frame(rs$rep.maov)
      dat$line_1 <- rownames(dat)
      dat <- melt(dat, id.vars = c("line_1"))
      colnames(dat) <- c("line_1", "line_2", "p_value")
      dat
    })    
    output$download_moav_all <- downloadHandler(
      filename = function() {"moav_results.csv"},
      content = function(file) {
        dat <- as.data.frame(rs$rep.maov)
        dat$line_1 <- rownames(dat)
        dat <- melt(dat, id.vars = c("line_1"))
        colnames(dat) <- c("line_1", "line_2", "p_value")
        write.csv(dat, file)
      }
    )    
    
    # AOV RESULTS 
    
    output$aov_results <- renderTable({
      if(is.null(rs$rep.aov)){return()}
      rep.aov <- rs$rep.aov
      rep.aov[(rep.aov$genotype_1 == input$to_plot & rep.aov$genotype_2 == input$ref_reps) | 
                       (rep.aov$genotype_1 == input$ref_reps & rep.aov$genotype_2 == input$to_plot),]
    })    
    output$download_oav <- downloadHandler(
      filename = function() {"oav_results.csv"},
      content = function(file) {
        rep.aov <- rs$rep.aov
        rep.aov[(rep.aov$genotype_1 == input$to_plot & rep.aov$genotype_2 == input$ref_reps) | 
                  (rep.aov$genotype_1 == input$ref_reps & rep.aov$genotype_2 == input$to_plot),]        
        write.csv(rep.aov, file)
      }
    )
    
    
    output$aov_results_all <- renderTable({
      if(is.null(rs$rep.aov)){return()}
      rs$rep.aov
    })    
    output$download_oav_all <- downloadHandler(
      filename = function() {"oav_results.csv"},
      content = function(file) {
        write.csv(rs$rep.aov, file)
      }
    )    
    
    
    output$aov_gene_results_all <- renderTable({
      if(is.null(rs$rep.aov.gene)){return()}
      rs$rep.aov.gene
    })    
    output$download_oav_gene_all <- downloadHandler(
      filename = function() {"oav_gene_results.csv"},
      content = function(file) {
        write.csv(rs$rep.aov.gene, file)
      }
    )    
    
    
    # GENE DATA
    
    output$gene_results_all <- renderTable({
      if(is.null(rs$gene)){return()}
      rs$gene
    })    
    output$download_gene_all <- downloadHandler(
      filename = function() {"gene_results.csv"},
      content = function(file) {
        write.csv(rs$gene, file)
      }
    )     
    
    
    # COMPARISON RESULTS 
    
    output$comp_results <- renderTable({
      if(is.null(rs$gene)){return()}
      rs$rep.mean
    })    
    output$download_comp <- downloadHandler(
      filename = function() {"matching_results.csv"},
      content = function(file) {
        write.csv(rs$rep.mean, file)
      }
    )    

})
