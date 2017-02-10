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
#   1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
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


##################################################################
#######   Plot functions
##################################################################


##### PLOT THE ROOT KEY IMAGE

plotRootKey <- function(root){
  
  #Initialize the root map
  root$tissue[root$tissue == "borders"] = NA
  
  pl <- ggplot() + coord_fixed() +
    theme_classic() +
    theme(
      axis.line = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size=18, face="bold"),
      #legend.position = "none",
      legend.text = element_text(size=18),
      legend.title = element_text(size=20),
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA)) 
  
  
  plot <- pl + geom_raster(data = root, aes(y, -x, fill = factor(tissue)), interpolate = T) + 
    scale_fill_brewer(palette = "Set1", na.value = "black")
  
  
  plot
}


##### PLOT THE ROOT REPORTER IMAGE

plotRootReporters <- function(reps, to_plot, root, rep.aov, 
                              rep.melt, rep.agg.short, sig, show){
  
  #Initialize the root map
  
  root$value <- NA
  root$value[root$tissue != "borders"] = 1
  
  n1 <- length(unique(rep.melt$root[rep.melt$line == reps]))
  n2 <- length(unique(rep.melt$root[rep.melt$line == to_plot]))
  
  # Get the value to plot for the reference line
  temp1 <- rep.agg.short[rep.agg.short$line == reps,]
  for(t in unique(as.character(temp1$variable))){
    root$value[root$id == 1 & root$tissue == t] <- temp1$value[temp1$variable == t]
  }
  
  # Get the value to plot for the comparison line
  temp2 <- rep.agg.short[rep.agg.short$line == to_plot ,]
  for(t in unique(as.character(temp2$variable))){
    root$value[root$id == 2 & root$tissue == t] <- temp2$value[temp1$variable == t]
  } 
  
  
  # Get the value to plot for the difference
  temp <- temp1
  temp$value <- temp2$value - temp1$value
  temp$pvalue <- rep.aov$pvalue[(rep.aov$genotype_1 == to_plot & rep.aov$genotype_2 == reps) | 
                                     (rep.aov$genotype_1 == reps & rep.aov$genotype_2 == to_plot)]
  temp <- temp[as.numeric(temp$pvalue) < 0.05,]
  
  root$value[root$id == 3 & root$tissue != "borders"] <- 0
  for(t in unique(as.character(temp$variable))){
    root$value[root$id == 3 & root$tissue == t] <- temp$value[temp$variable == t]
  }             
  
  root1 <- root[root$id == 1,]
  root2 <- root[root$id == 2,]
  root3 <- root[root$id == 3,]
  
  rg <- range(root1$value, root2$value, na.rm = T)
  
  pl <- ggplot() + coord_fixed() +
    theme_classic() +
    theme(
      axis.line = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size=12, face="bold"),
      #legend.position = "none",
      legend.text = element_text(size=14),
      legend.title = element_text(size=15),
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA)) 
  
  
  plot1 <- pl + geom_raster(data = root1, aes(y, -x, fill = value), interpolate = T) + 
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=c("#FFFFFF", "#51A9F9","#255A8A","#000000"),na.value = "black", limits=rg) + 
    ggtitle(paste0(reps, "\n n = ",n1))
  
  plot2 <- pl + geom_raster(data = root2, aes(y, -x, fill = value), interpolate = T) + 
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=c("#FFFFFF", "#51A9F9","#255A8A","#000000"),na.value = "black", limits=rg) +     
    ggtitle(paste0(to_plot, "\n n = ",n2)) 
  
  if(show){
    lim <- c(-max(abs(range(root3$value, na.rm = T))), max(abs(range(root3$value, na.rm = T))))
    
    ti <- "SIGNIFICANT\nDIFFERENCES"
    if(sig) ti <- "SIGNIFICANT\nDIFFERENCES\n***"
    plot3 <- pl + geom_raster(data = root3, aes(y, -x, fill = value), interpolate = T) + 
      scale_fill_gradientn(colours=c("#EF7C09", "#FFFFFF","#51A9F9"), limits = lim, na.value = "black") + 
      ggtitle(ti)
    
    plot <- grid.arrange(plot1, plot2, plot3, ncol=3)
  }else{
    plot <- grid.arrange(plot1, plot2, ncol=2)
  }
  
  plot
}





##### PLOT THE ROOT GENE IMAGE

plotRootGene <- function(reps, gene, root, rep.agg.short){
  
  #Initialize the root map
  
  root$value <- NA
  root$value[root$tissue != "borders"] = 1

  # Get the value to plot for the reference line
  temp1 <- rep.agg.short[rep.agg.short$line == reps,]
  name <- temp1$line[1]
  match <- temp1$match[1]
  for(t in unique(as.character(temp1$variable))){
    root$value[root$id == 1 & root$tissue == t] <- temp1$value[temp1$variable == t]
  }
  
  # Get the value to plot for the comparison line
  temp <- gene[gene$Gene_ID == match,]
  temp <- melt(temp, id.vars = c("Gene_ID", "match", "distance"))  
  for(t in unique(as.character(temp$variable))){
    root$value[root$id == 2 & root$tissue == t] <- temp$value[temp$variable == t]
  }   
    
  root1 <- root[root$id == 1,]
  root2 <- root[root$id == 2,]

  rg <- range(root1$value, root2$value, na.rm = T)
  
  pl <- ggplot() + coord_fixed() +
    theme_classic() +
    theme(
      axis.line = element_blank(), 
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size=12, face="bold"),
      #legend.position = "none",
      legend.text = element_text(size=14),
      legend.title = element_text(size=15),
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA)) 
  
  
  plot1 <- pl + geom_raster(data = root1, aes(y, -x, fill = value), interpolate = T) + 
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=c("#FFFFFF", "#51A9F9","#255A8A","#000000"),na.value = "black", limits=rg) +     
    ggtitle(name)
  
  plot2 <- pl + geom_raster(data = root2, aes(y, -x, fill = value), interpolate = T) + 
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=c("#FFFFFF", "#51A9F9","#255A8A","#000000"),na.value = "black", limits=rg) +     
    ggtitle(match) 
  
  plot <- grid.arrange(plot1, plot2, ncol=2)
  
  plot
}




###### BARPLOTS

barplot_comp_1 <- function(reps, gene, rep.agg.short){
  
  temp2 <- rep.agg.short[rep.agg.short$line == reps,]
  match <- temp2$match[1]
  temp <- melt(gene[gene$Gene_ID == match,], id.vars = c("Gene_ID", "match", "distance"))
  temp <- temp[,-2]
  colnames(temp) <- colnames(temp2)
  temp <- rbind(temp, temp2)
  
  plot1 <- ggplot(temp, aes(variable, value, fill=line)) + 
    geom_bar(stat = "identity", position=position_dodge(width=0.9), width=0.8) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=14),
      axis.text = element_text(size=20),
      axis.title =  element_text(size=20)
    ) +
    xlab("") + ylab("Relative expression value")
  #p <- ggplotly(plot1)
  plot1
  
}


barplot_comp <- function(reps, to_plot, data){

  temp <- data[data$line == reps | data$line == to_plot,]
  
  plot1 <- ggplot(temp, aes(variable, value, colour=line)) + 
    geom_boxplot( width=0.8, size=1.5) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=14),
      axis.text = element_text(size=20),
      axis.title =  element_text(size=20)
    ) +
    xlab("") + ylab("Relative expression value")

    plot1
  
}




## HEATMAP #############################

heatmap <- function(rep.maov){
  
  rep.maov[rep.maov == 0] <- "ns."
  rep.maov[rep.maov <= 0.01] <- "***"
  rep.maov[rep.maov <= 0.05 & rep.maov > 0.01] <- "*"
  rep.maov[rep.maov > 0.05] <- "ns."
  
  dat <- as.data.frame(rep.maov)
  dat$line_1 <- rownames(dat)
  dat <- melt(dat, id.vars = c("line_1"))
  
  dat$line_2 <- dat$variable
  dat$p_value <- dat$value
  ## Example data
  plot1 <- ggplot(dat, aes(line_1, line_2, z= p_value)) + 
    geom_tile(aes(fill = p_value)) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=15),
      axis.text = element_text(size=15),
      legend.text = element_text(size=30),
      legend.title = element_text(size=30)) +
    
    scale_fill_manual(values=c("#00BC4760", "#00BC47", "white"), 
                      name="Significance level",
                      labels=c("*","***", "ns.")) +
    #scale_fill_gradient(low="blue", high="white", space="Lab")   +
    xlab("") + ylab("") + 
    coord_fixed()
  
  plot1
}


heatmap_dist <- function(gene.dist){
  
  dat <- as.data.frame(t(gene.dist))
  dat$line_1 <- rownames(dat)
  dat <- melt(dat, id.vars = c("line_1"))
  
  dat$line_2 <- dat$variable
  dat$distance <- dat$value
  ## Example data
  plot1 <- ggplot(dat, aes(line_1, line_2, z= distance)) + 
    geom_tile(aes(fill = distance)) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=15),
      axis.text = element_text(size=15),
      legend.text = element_text(size=30),
      legend.title = element_text(size=30)) +
    scale_fill_distiller(palette="Spectral")+
    #scale_fill_gradient(low="blue", high="white", space="Lab")   +
    xlab("") + ylab("")  +
    coord_fixed()
  
  plot1
}


## LDA PLOT #############################

ldaPlot <- function(reps, to_plot, rep){
  
  rep$col <- 0
  rep$alpha <- 0.2
  rep$size <- 1
  rep$col[rep$line == to_plot] <- 1
  rep$col[rep$line == reps] <- 2
  rep$alpha[rep$line == to_plot] <- 1
  rep$alpha[rep$line == reps] <- 1
  
  
  plot2 <- ggplot() + 
    geom_point(data=rep[rep$col==0,], aes(PC1, PC2, fill=line, alpha=alpha), shape=1) +
    scale_fill_grey() + 
    geom_point(data=rep[rep$col>0,], aes(PC1, PC2, colour=line), shape=16) + 
    stat_ellipse(data=rep[rep$col>0,], aes(PC1, PC2, colour=line), level = 0.6, size=1) + 
    theme_bw() + 
    theme(legend.position="none") + 
    scale_colour_manual(values=c("#F96F70", "#00BC47"), 
                        name="Lines",
                        labels=c(to_plot, reps)) 
  
  
  p <- ggplotly(plot2)
  p
}
