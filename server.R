
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

    # load("www/example/root.RData")
    root <- fread("www/example/root.csv")
    remove <- c(65,146,116,416,443,454,490)
    root <- root[!(root$id %in% remove),]
    rs <- reactiveValues(reporter = NULL, 
                         gene = NULL,
                         lda.fit = NULL,
                         rep.melt = NULL,
                         tissues = NULL,
                         rep.aov = NULL,
                         rep.aov.gene = NULL,
                         rep.fit.gene = NULL,
                         rep.maov = NULL,
                         rep.fit = NULL,
                         rep.fit.table = NULL,
                         rep.pearson = NULL,
                         rep.spearman = NULL,
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

      # fileName <- 'www/datasets/reporters/'
      # flist <- unique(gsub(".csv", "", list.files(fileName)))
      flist <-c("experession_full", "experession_ji_young")
      fl <- dts$name[dts$id %in% flist]
      updateSelectInput(session, "reporters", choices = fl)  
      
      # fileName <- 'www/datasets/microarrays/'
      # flist <- unique(gsub(".csv", "", list.files(fileName)))
      flist <-c("all_genes")
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

         # gene <- fread("www/datasets/microarrays/all_genes.csv")
         # temp <- read_rsml("www/datasets/reporters/experession_ji_young.rsml")
         # temp <- read_rsml("~/Downloads/reporter.rsml")
         # temp <- read_rsml("www/datasets/reporters/experession_full.rsml")
         
          # Load two datafiles, for the gene and the reporters

          if(input$use_example){
            dat1 <- rs$dts$id[rs$dts$name == input$microarrays]
            dat2 <- rs$dts$id[rs$dts$name == input$reporters]
            if(dat1 == "all_genes" & dat2 == "experession_full"){
              if(input$log2bis){
                message("Loading dataset 1 - log")
                load("www/example/example_log1.RData")
              }else{
                message("Loading dataset 1")
                load("www/example/example1.RData")
              }
            }else if(dat1 == "all_genes" & dat2 == "experession_ji_young"){
              if(input$log2bis){
                message("Loading dataset 2 - log")
                load("www/example/example_log2.RData")
              }else{
                message("Loading dataset 2")
                load("www/example/example2.RData")
              }
            }            
            
          }else{
            inGene <- input$gene_file
            inReporter <- input$rep_file

            if (is.null(inReporter)) return(NULL)
            if(!is.null(inGene)){
              gene <- fread(inGene$datapath)[,-4]   
            }
            temp <- read_rsml(inReporter$datapath)   
          
            temp <- temp[temp$cell_type %in% input$type_to_analyse,]

            # Transform the data
            if(input$log2) temp$value <- log2(temp$value)
            if(!input$use_absolute) temp$value <- ddply(temp, .(line, root), plyr::summarize, value=scale(value))$value
            
            # Average the data by line, root, cell type
            if(input$method == "Mean") mean_data <- ddply(temp, .(line, root, cell_type), plyr::summarize, value=mean(value))
            if(input$method == "Median") mean_data <- ddply(temp, .(line, root, cell_type), plyr::summarize, value=median(value))
            if(input$method == "Min") mean_data <- ddply(temp, .(line, root, cell_type), plyr::summarize, value=min(value))
            if(input$method == "Max") mean_data <- ddply(temp, .(line, root, cell_type), plyr::summarize, value=max(value))
            mean_data <- mean_data[!is.na(mean_data$value),]
            mean_data <- mean_data[!is.infinite(mean_data$value),]
            
            # Average the data by line, cell type
            mean_data_2 <- ddply(mean_data, .(line, cell_type), plyr::summarise, value=mean(value))
            mean_data_2 <- mean_data_2[!is.na(mean_data_2$value),]          
            mean_data_2 <- mean_data_2[!is.infinite(mean_data_2$value),]          
              
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
            rep.agg.short <- ddply(rep.melt, .(line, variable), plyr::summarize, value=mean(value))
            
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
            
            fit.results <<- matrix(NA, ncol = l1, nrow = l1) 
            colnames(fit.results) <- p.list 
            rownames(fit.results) <- p.list
            
            pearson.results <<- matrix(NA, ncol = l1, nrow = l1) 
            colnames(pearson.results) <- p.list 
            rownames(pearson.results) <- p.list
            
            spearman.results <<- matrix(NA, ncol = l1, nrow = l1) 
            colnames(spearman.results) <- p.list 
            rownames(spearman.results) <- p.list
              
            fit.table <- NULL
            
            i <- 1
            k <- 1
            for(p in p.list){
              j <- 1
              for(p1 in p.list){
                if(p != p1 & is.na(maov.results[j,i])){
                  
                  # MANOVA analysis to compare the plants
                  temp <- reporter[reporter$line == p | reporter$line == p1,]
                  
                  # Get the ines sleected by the user
                  # ts <- c("columella","cortex","endodermis","epidermis","lateralrootcap","QC","stele")
                  ts <- input$type_to_analyse
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
  
                  
                  temp2 <- melt(temp, id=c("line", "root"))
                  temp2 <- ddply(temp2, .(line, variable), plyr::summarise, value=mean(value))
                  x <- temp2$value[temp2$line == p]; y <- temp2$value[temp2$line == p1] 
                  
                  fit.results[j,i] <- summary(lm(y ~ x))$r.squared
                  fit.results[i,j] <- fit.results[j,i]
                  
                  pearson.results[j,i] <- rcorr(x, y, type = "pearson")[[1]][1,2]
                  pearson.results[i,j] <- pearson.results[j,i]
                  
                  spearman.results[j,i] <- rcorr(x, y, type = "spearman")[[1]][1,2]
                  spearman.results[i,j] <- spearman.results[j,i]
                  
                  
                  fit.table <- rbind(fit.table, data.table(line1 = p, line2=p1, 
                                                           r2 = fit.results[j,i], 
                                                           r.pvalue = summary(lm(y ~ x))$coefficients[2, 4],
                                                           pearson = pearson.results[j,i],
                                                           pearson.pvalue = rcorr(x, y, type = "pearson")[[3]][1,2],
                                                           spearman = spearman.results[j,i],
                                                           spearman.pvalue = rcorr(x, y, type = "spearman")[[3]][1,2]))
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
          # ----- AOV BETWEEN THE REPORTER AND THE CORRESPONDING GENE   ---------------------------
          # ----------------------------------------------------------------------------------            
          
            incProgress(1/5, detail = "Matching the genes")
            
            aov.results.gene <- NULL
            fit.results.gene <- NULL
            
            for(p in p.list){
              ge <- strsplit(p, "_")[[1]]
              ge <- ge[grepl("AT[1-9]", ge)] # find the corresping gene
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
                
                temp <- ddply(temp, .(line, variable), plyr::summarise, value=mean(value))
                temp2 <- ddply(temp2, .(line, variable), plyr::summarise, value=mean(value))
                temp <- merge(temp, temp2, by.x="variable", by.y="variable")
                fit.result <- summary(lm(temp$value.y ~ temp$value.x))$r.squared
                pearson.result <- rcorr(temp$value.x, temp$value.y, type = "pearson")[[1]][1,2]
                spearman.result <- rcorr(temp$value.x, temp$value.y, type = "spearman")[[1]][1,2]    
                
                fit.results.gene <- rbind(fit.results.gene, data.frame(line=p, gene=ge, 
                                                                       r2=fit.result, r.pvalue = summary(lm(temp$value.y ~ temp$value.x))$coefficients[2, 4],
                                                                       pearson=pearson.result, pearson.pvalue = rcorr(temp$value.x, temp$value.y, type = "pearson")[[3]][1,2], 
                                                                       spearman = spearman.result, spearman.pvalue = rcorr(temp$value.x, temp$value.y, type = "spearman")[[3]][1,2]))                
              }
            }
            
          
          
            # -------------------------------------------------------
            # ----- PCA ANALYSIS ------------------------------------
            # -------------------------------------------------------            
          
            incProgress(1/5, detail = "Performing the PCA")
          
            # Make the LDA analysis on the reporter dataset.
            # ts <- c("line", input$type_to_analyse)
            temp <- reporter[complete.cases(reporter),]
            lines <- temp$line
            temp <- temp[,input$type_to_analyse]
            # temp <- temp[,cell_types]
            
            
            
            pca <- prcomp(temp, retx = T, scale=T)  # Make the PCA
            pca.results <- cbind(line=lines, data.frame(pca$x)[,])

          }
        
          # save(reporter,rep.mean,gene,pca.results,rep.melt,tissues,aov.results,aov.results.gene,maov.results,p.list,rep.agg.short, file="www/example/example2.RData")
        
          # rs$reporter <- reporter
          rs$rep.mean <- rep.mean
          rs$gene <- gene
          #rs$lda.fit <- fit
          rs$lda.fit <- pca.results
          rs$reporter <- reporter
          rs$rep.melt <- rep.melt
          rs$tissues <- tissues
          rs$rep.aov <- aov.results
          rs$rep.aov.gene <- aov.results.gene
          rs$rep.fit.gene <- fit.results.gene
          rs$rep.maov <- maov.results
          rs$rep.fit <- fit.results
          rs$rep.fit.table <- fit.table
          rs$rep.pearson <- pearson.results
          rs$rep.spearman <- spearman.results
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
      
      if(sel %in% reps){
        message(">>>>>>>>>>>>")
        message(sel)
        message(">>>>>>>>>>>>")
      }else{ sel <- reps1[1]}

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
      ct_options <- list()
      sel <- input$type_to_plot
      if(length(sel) == 0) sel = cell_types
      if(!is.null(input$type_to_analyse)) sel <- sel[sel %in% input$type_to_analyse]
      for(ct in cell_types) ct_options[[ct]] <- ct
      updateSelectInput(session, "type_to_plot", choices = ct_options, selected=sel) 
    })  
    
    observe({
      ct_options <- list()
      sel <- input$type_to_plot_2
      if(length(sel) == 0) sel = cell_types
      if(!is.null(input$type_to_analyse)) sel <- sel[sel %in% input$type_to_analyse]
      for(ct in cell_types) ct_options[[ct]] <- ct
      updateSelectInput(session, "type_to_plot_2", choices = ct_options, selected=sel) 
    })  
    
    observe({
      
      if(is.null(rs$reporter) || grepl("Load", input$ref_reps)){return()}
      
      sel <- input$ref_reps
      reps <- na.omit(unique(factor(rs$reporter$line)))
      reps2 <- reps[-match(input$to_plot, reps)]
      
      if(sel %in% reps){
      }else{ sel <- reps2[1]}
      
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

    
## FITPLOT #############################
    
    output$fitPlot <- renderPlot({
      if(is.null(rs$reporter)){return()}
      
      temp <- rs$reporter[rs$reporter$line == input$ref_reps | rs$reporter$line == input$to_plot,]
      temp2 <- melt(temp, id=c("line", "root"))
      
      print(fitPlot(dat = temp2,
                           p = input$ref_reps, 
                           p1 = input$to_plot,
                          types = input$type_to_plot
      ))
    })
    
## FITPLOT GENE #############################
    
    output$fitPlot_1 <- renderPlot({
      if(is.null(rs$gene)){return()}
      
      temp2 <- rs$rep.melt[rs$rep.melt$line == input$ref_reps_2,]
      match <- temp2$match[1]
      temp <- rs$gene[rs$gene$Gene_ID == match,]
      # temp2 <- temp2[,-c(2, 5)]
      # colnames(temp) <- colnames(temp2)
      # temp <- rbind(temp, temp2)
      
      print(fitPlot_1(dat = temp, dat1 = temp2,
                    p = input$ref_reps_2, 
                    p1 = match,
                    types = input$type_to_plot_2
      ))
    })    
    
        

## BARPLOT #############################
    
    output$barplot_comp_1 <- renderPlot({
      if(is.null(rs$gene)){return()}
      print(barplot_comp_1(reps = input$ref_reps_2,
                           rep.agg = rs$rep.melt, 
                           gene = rs$gene,
                           types = input$type_to_plot_2
                            ))
    })  
 
    
## BARPLOT #############################
    
    output$barplot_comp <- renderPlot({
      if(is.null(rs$reporter)){return()}
      print(barplot_comp(input$ref_reps, input$to_plot, rs$rep.melt, types = input$type_to_plot))
    })   

## PLOT THE ROOT #############################
    
    output$plotRootKey <- renderPlot({
      print(plotRootKey(rs$root))
    }) 

    
## PLOT THE ROOT WITH THE REPORTER DATA  #############################
    
    output$plotRoot <- renderPlot({
      if(is.null(rs$reporter)){return()}
      abs = F
      if(!is.null(input$use_absolute)) abs <- input$use_absolute
      print(plotRootReporters(reps = input$ref_reps, 
                              to_plot = input$to_plot, 
                              rep.aov = rs$rep.aov, 
                              root = rs$root, 
                              rep.melt = rs$rep.melt, 
                              rep.agg.short = rs$rep.agg.short,
                              sig = rs$rep.maov[input$ref_reps, input$to_plot] < 0.05,
                              input$show_diff,
                              range = c(0,1),#input$display_range,
                              types = input$type_to_plot,
                              abs = abs
                              ))
    }) 
    
## PLOT THE ROOT WITH THE GENE DATA  #############################
    
    output$plotRootGene <- renderPlot({
      if(is.null(rs$gene)){return()}
      abs = F
      if(!is.null(input$use_absolute)) abs <- input$use_absolute
      print(plotRootGene(reps = input$ref_reps_2, 
                         root = rs$root, 
                         gene = rs$gene,
                         rep.agg.short = rs$rep.agg.short,
                         range = c(0,1),#input$display_range_1,
                         types = input$type_to_plot_2,
                         abs = abs
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
    
    
## PLOT THE HEATMAP OF THE REGRESSIONS ANALYSIS  
    
    output$heatmap_fit <- renderPlot({
      if(is.null(rs$rep.maov)){return()}
      print(heatmap_fit(rs$rep.fit))
    }) 
    
    
## PLOT THE HEATMAP OF THE PEARSON ANALYSIS  
    
    output$heatmap_spearman <- renderPlot({
      if(is.null(rs$rep.spearman)){return()}
      print(heatmap_fit(rs$rep.spearman, diff=T))
    })     
    
## PLOT THE HEATMAP OF THE PEARSON ANALYSIS  
    
    output$heatmap_pearson <- renderPlot({
      if(is.null(rs$rep.pearson)){return()}
      print(heatmap_fit(rs$rep.pearson, diff=T))
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
    
    
    
    output$line_comp_fit <- renderText({ 
      if(is.null(rs$rep.fit) || grepl("Load", input$to_plot)){return()}
      
      sig <- rs$rep.fit[input$ref_reps, input$to_plot]
      
      # Update the r/squared value based on the selected tissue types
      temp <- rs$reporter[rs$reporter$line == input$ref_reps | rs$reporter$line == input$to_plot,]
      temp2 <- melt(temp, id=c("line", "root"))
      temp2 <- ddply(temp2, .(line, variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp2 <- temp2[temp2$variable %in% input$type_to_plot,]
      x <- temp2$avg[temp2$line == input$ref_reps]; 
      y <- temp2$avg[temp2$line == input$to_plot] 
      sig.new <- summary(lm(y~x))$r.squared
      
      text <- ""
      if(!is.na(sig)){
        text <- paste0("r-squared: All = ",round(sig, 4), " || Selected = ",round(sig.new, 4))
      }
      return(HTML(text))
    })   
    
    output$line_comp_pears <- renderText({ 
      if(is.null(rs$rep.pearson) || grepl("Load", input$to_plot)){return()}
      sig <- rs$rep.pearson[input$ref_reps, input$to_plot]
      
      # Update the r/squared value based on the selected tissue types
      temp <- rs$reporter[rs$reporter$line == input$ref_reps | rs$reporter$line == input$to_plot,]
      temp2 <- melt(temp, id=c("line", "root"))
      temp2 <- ddply(temp2, .(line, variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp2 <- temp2[temp2$variable %in% input$type_to_plot,]
      x <- temp2$avg[temp2$line == input$ref_reps]; 
      y <- temp2$avg[temp2$line == input$to_plot] 
      sig.new <- rcorr(x, y, type = "pearson")[[1]][1,2]
      
      text <- ""
      if(!is.na(sig)){
        text <- paste0("Pearson correlation: All = ",round(sig, 4), " || Selected = ",round(sig.new, 4))
      }
      return(HTML(text))
    })   
    
    output$line_comp_spear <- renderText({ 
      if(is.null(rs$rep.spearman) || grepl("Load", input$to_plot)){return()}
      sig <- rs$rep.spearman[input$ref_reps, input$to_plot]
      
      # Update the r/squared value based on the selected tissue types
      temp <- rs$reporter[rs$reporter$line == input$ref_reps | rs$reporter$line == input$to_plot,]
      temp2 <- melt(temp, id=c("line", "root"))
      temp2 <- ddply(temp2, .(line, variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp2 <- temp2[temp2$variable %in% input$type_to_plot,]
      x <- temp2$avg[temp2$line == input$ref_reps]; 
      y <- temp2$avg[temp2$line == input$to_plot] 
      sig.new <- rcorr(x, y, type = "spearman")[[1]][1,2]
      
      text <- ""
      if(!is.na(sig)){
        text <- paste0("Spearman correlation: All = ",round(sig, 4), " || Selected = ",round(sig.new, 4))
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
    
    output$gene_comp_fit <- renderText({ 
      if(is.null(rs$rep.fit.gene) || grepl("Load", input$to_plot)){return()}
      text <- ""
      sig <- rs$rep.fit.gene[rs$rep.fit.gene$line == input$ref_reps_2, c("r2")]
      
      temp2 <- rs$rep.melt[rs$rep.melt$line == input$ref_reps_2,]
      match <- temp2$match[1]
      temp <- rs$gene[rs$gene$Gene_ID == match,]
      temp1 <- ddply(temp, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp2 <- ddply(temp2, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp <- merge(temp1, temp2, by="variable")
      temp <- temp[temp$variable %in% input$type_to_plot_2,]
      x <- temp$avg.x 
      y <- temp$avg.y
      sig.new <- summary(lm(y~x))$r.squared

      if(!is.na(sig)){
        text <- paste0("R-squared: All = ",round(sig, 4), " || Selected = ",round(sig.new, 4))
      }
      return(HTML(text))
    }) 
    
    output$gene_comp_pears <- renderText({ 
      if(is.null(rs$rep.fit.gene) || grepl("Load", input$to_plot)){return()}
      text <- ""
      sig <- rs$rep.fit.gene[rs$rep.fit.gene$line == input$ref_reps_2, c("pearson")]
      
      temp2 <- rs$rep.melt[rs$rep.melt$line == input$ref_reps_2,]
      match <- temp2$match[1]
      temp <- rs$gene[rs$gene$Gene_ID == match,]
      temp1 <- ddply(temp, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp2 <- ddply(temp2, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp <- merge(temp1, temp2, by="variable")
      temp <- temp[temp$variable %in% input$type_to_plot_2,]
      x <- temp$avg.x 
      y <- temp$avg.y
      sig.new <- rcorr(x, y, type = "pearson")[[1]][1,2]
      
      if(!is.na(sig)){
        text <- paste0("Pearson correlation: All = ",round(sig, 4), " || Selected = ",round(sig.new, 4))
      }
      return(HTML(text))
    }) 
    
    output$gene_comp_spear <- renderText({ 
      if(is.null(rs$rep.fit.gene) || grepl("Load", input$to_plot)){return()}
      text <- ""
      sig <- rs$rep.fit.gene[rs$rep.fit.gene$line == input$ref_reps_2, c("spearman")]
      
      temp2 <- rs$rep.melt[rs$rep.melt$line == input$ref_reps_2,]
      match <- temp2$match[1]
      temp <- rs$gene[rs$gene$Gene_ID == match,]
      temp1 <- ddply(temp, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp2 <- ddply(temp2, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
      temp <- merge(temp1, temp2, by="variable")
      temp <- temp[temp$variable %in% input$type_to_plot_2,]
      x <- temp$avg.x 
      y <- temp$avg.y
      sig.new <- rcorr(x, y, type = "spearman")[[1]][1,2]
      
      
      if(!is.na(sig)){
        text <- paste0("Spearman correlation: All = ",round(sig, 4), " || Selected = ",round(sig.new, 4))
      }
      return(HTML(text))
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
    

    output$fit_results_all <- renderTable({
      if(is.null(rs$rep.fit)){return()}
      dat <- as.data.frame(rs$rep.fit.table)
      # dat$line_1 <- rownames(dat)
      # dat <- melt(dat, id.vars = c("line_1"))
      # colnames(dat) <- c("line_1", "line_2", "r2")
      dat
    })    
    output$download_fit_all <- downloadHandler(
      filename = function() {"moav_results.csv"},
      content = function(file) {
        dat <- as.data.frame(rs$rep.fit.table)
        # dat$line_1 <- rownames(dat)
        # dat <- melt(dat, id.vars = c("line_1"))
        # colnames(dat) <- c("line_1", "line_2", "r2")
        write.csv(dat, file)
      }
    )
    
    output$spear_results_all <- renderTable({
      if(is.null(rs$rep.spearman)){return()}
      dat <- as.data.frame(rs$rep.spearman)
      dat$line_1 <- rownames(dat)
      dat <- melt(dat, id.vars = c("line_1"))
      colnames(dat) <- c("line_1", "line_2", "correlation")
      dat
    })    
    output$download_spear_all <- downloadHandler(
      filename = function() {"spearman_results.csv"},
      content = function(file) {
        dat <- as.data.frame(rs$rep.spearman)
        dat$line_1 <- rownames(dat)
        dat <- melt(dat, id.vars = c("line_1"))
        colnames(dat) <- c("line_1", "line_2", "correlation")
        write.csv(dat, file)
      }
    )
    
    
    output$pears_results_all <- renderTable({
      if(is.null(rs$rep.pearson)){return()}
      dat <- as.data.frame(rs$rep.maov)
      dat$line_1 <- rownames(dat)
      dat <- melt(dat, id.vars = c("line_1"))
      colnames(dat) <- c("line_1", "line_2", "correlation")
      dat
    })    
    output$download_pears_all <- downloadHandler(
      filename = function() {"pearson_results.csv"},
      content = function(file) {
        dat <- as.data.frame(rs$rep.pearson)
        dat$line_1 <- rownames(dat)
        dat <- melt(dat, id.vars = c("line_1"))
        colnames(dat) <- c("line_1", "line_2", "correlation")
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
    
    output$fit_gene_results_all <- renderTable({
      if(is.null(rs$rep.aov.gene)){return()}
      rs$rep.fit.gene
    })    
    output$download_fit_gene_all <- downloadHandler(
      filename = function() {"fit_gene_results.csv"},
      content = function(file) {
        write.csv(rs$rep.fit.gene, file)
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
