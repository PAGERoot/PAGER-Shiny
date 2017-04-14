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
  pl <- ggplot(root, aes(x=x, y=y)) +
    geom_polygon(aes(fill=tissue, group=id), colour="black") +
    coord_fixed()+
    theme_classic() +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5, size=18, face="bold"))
  
  pl
}


##### PLOT THE ROOT REPORTER IMAGE

plotRootReporters <- function(reps, to_plot, root, rep.aov, 
                              rep.melt, rep.agg.short, sig, 
                              show, range, types){
  
  #Initialize the root map
  root$value <- "NA"
  
  root1 <- root
  root2 <- root
  root3 <- root
  
  n1 <- length(unique(rep.melt$root[rep.melt$line == reps]))
  n2 <- length(unique(rep.melt$root[rep.melt$line == to_plot]))
  
  # Get the value to plot for the reference line
  temp1 <- rep.agg.short[rep.agg.short$line == reps,]
  temp1 <- temp1[temp1$variable %in% types,]
  temp1$value <- range01(temp1$value)
  for(t in unique(as.character(temp1$variable))){
    root1$value[root1$tissue == t] <- temp1$value[temp1$variable == t]
  }
  
  
  
  # Get the value to plot for the comparison line
  temp2 <- rep.agg.short[rep.agg.short$line == to_plot ,]
  temp2 <- temp2[temp2$variable %in% types,]
  temp2$value <- range01(temp2$value)
  for(t in unique(as.character(temp2$variable))){
    root2$value[root2$tissue == t] <- temp2$value[temp1$variable == t]
  } 
  
  
  
  # Get the value to plot for the difference
  temp <- temp1
  for(t in unique(as.character(temp2$variable))){
    temp$value[temp$variable == t] <- temp2$value[temp2$variable == t] - temp1$value[temp1$variable == t]
  }
  pvals <- rep.aov[(rep.aov$genotype_1 == to_plot & rep.aov$genotype_2 == reps) | 
                                     (rep.aov$genotype_1 == reps & rep.aov$genotype_2 == to_plot),c("tissue", "pvalue")]
  temp <- merge(temp, pvals, by.x="variable", by.y = "tissue")  
  temp <- temp[as.numeric(temp$pvalue) < 0.05,]
  
  for(t in unique(as.character(temp$variable))){
    root3$value[root3$tissue == t] <- temp$value[temp$variable == t]
  }             
  
  rg <- range
  # rg <- range(0,1)
  
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
      legend.text = element_text(size=12),
      legend.title = element_text(size=12),
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA))  + 
    geom_rect(aes(ymin=-200, ymax=-300, xmin=min(root$x), xmax=max(root$x)+10), fill="white")
  
  
  plot1 <- pl + geom_polygon(data = root1, aes(x=x, y=y, fill=value, group=id), colour="black") +
    # geom_raster(data = root1, aes(y, -x, fill = value), interpolate = T) + 
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=cscale,na.value = "grey", limits=rg) + 
    ggtitle(paste0(reps, "\n n = ",n1))

  plot2 <- pl + geom_polygon(data = root2, aes(x=x, y=y, fill=value, group=id), colour="black") +
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours= cscale,na.value = "grey", limits=rg) +     
    ggtitle(paste0(to_plot, "\n n = ",n2))   

  if(show){
    lim <- c(-max(abs(range(root3$value, na.rm = T))), max(abs(range(root3$value, na.rm = T))))
    
    ti <- "SIGNIFICANT\nDIFFERENCES"
    if(sig) ti <- "SIGNIFICANT\nDIFFERENCES\n***"
    plot3 <- pl + geom_polygon(data = root3, aes(x=x, y=y, fill=value, group=id), colour="black") +
      scale_fill_gradientn(colours=c("#EF7C09", "#FFFFFF","#51A9F9"), limits = lim, na.value = "grey") + 
      ggtitle(ti)  

    plot <- grid.arrange(plot1, plot2, plot3, ncol=3)
  }else{
    plot <- grid.arrange(plot1, plot2, ncol=2)
  }
  
  plot
}





##### PLOT THE ROOT GENE IMAGE

plotRootGene <- function(reps, gene, root, rep.agg.short, range, types){
  
  #Initialize the root map
  root$value <- NA
  root1 <- root
  root2 <- root

  # Get the value to plot for the reference line
  temp1 <- rep.agg.short[rep.agg.short$line == reps,]
  temp1 <- temp1[temp1$variable %in% types,]
  temp1$value <- range01(temp1$value)
  name <- temp1$line[1]
  match <- temp1$match[1]
  for(t in unique(as.character(temp1$variable))){
    root1$value[root1$tissue == t] <- temp1$value[temp1$variable == t]
  }
  
  # Get the value to plot for the comparison line
  temp <- gene[gene$Gene_ID == match,]
  temp <- ddply(temp, .(Gene_ID, variable), summarise, value=mean(value))#melt(temp, id.vars = c("Gene_ID"))  
  temp <- temp[temp$variable %in% types,]
  temp$value <- range01(temp$value)
  for(t in unique(as.character(temp$variable))){
    root2$value[root2$tissue == t] <- temp$value[temp$variable == t]
  }   

  rg <- range
  
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
      legend.text = element_text(size=12),
      legend.title = element_text(size=12),
      panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
      plot.background = element_rect(fill = "transparent",colour = NA)) 
  
  plot1 <- pl + geom_polygon(data = root1, aes(x=x, y=y, fill=value, group=id), colour="black") +
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=cscale,na.value = "grey", limits=rg) +     
    ggtitle(name)  
  
  plot2 <- pl + geom_polygon(data = root2, aes(x=x, y=y, fill=value, group=id), colour="black") +
    # scale_fill_gradientn(colors = terrain.colors(7), na.value = "black", limits=rg) + 
    scale_fill_gradientn(colours=cscale,na.value = "grey", limits=rg) +     
    ggtitle(match)  
  
  plot <- grid.arrange(plot1, plot2, ncol=2)
  
  plot
}




###### BARPLOTS

barplot_comp_1 <- function(reps, gene, rep.agg, types){
  
  temp2 <-rep.agg[rep.agg$line == reps,]
  match <- temp2$match[1]
  temp <- gene[gene$Gene_ID == match,]
  temp2 <- temp2[,-c(2, 5)]
  colnames(temp) <- colnames(temp2)
  temp <- rbind(temp, temp2)
  
  temp <- temp[temp$variable %in% types,]
  
  temp$variable <- as.character(temp$variable)
  temp$variable[temp$variable == "columella"] <- "01 - Columella"
  temp$variable[temp$variable == "lateralrootcap"] <- "02 - Lateral root cap"
  temp$variable[temp$variable == "QC"] <- "03 - Quiecent center"
  temp$variable[temp$variable == "epidermis"] <- "04 - Epidermis"
  temp$variable[temp$variable == "cortex"] <- "05 - Cortex"
  temp$variable[temp$variable == "endodermis"] <- "06 - Endodermis"
  temp$variable[temp$variable == "stele"] <- "07 - Stele"  
  
  plot1 <- ggplot(temp, aes(factor(variable), value, fill=line)) + 
    # geom_bar(stat = "identity", position=position_dodge(width=0.9), width=0.8) + 
    geom_boxplot( width=0.8, size=1) + 
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


barplot_comp <- function(reps, to_plot, data, types){

  temp <- data[data$line == reps | data$line == to_plot,]
  
  temp <- temp[temp$variable %in% types,]
  
  
  temp$variable <- as.character(temp$variable)
  temp$variable[temp$variable == "columella"] <- "01 - Columella"
  temp$variable[temp$variable == "lateralrootcap"] <- "02 - Lateral root cap"
  temp$variable[temp$variable == "QC"] <- "03 - Quiecent center"
  temp$variable[temp$variable == "epidermis"] <- "04 - Epidermis"
  temp$variable[temp$variable == "cortex"] <- "05 - Cortex"
  temp$variable[temp$variable == "endodermis"] <- "06 - Endodermis"
  temp$variable[temp$variable == "stele"] <- "07 - Stele"
  
  plot1 <- ggplot(temp, aes(factor(variable), value, fill=line)) + 
    geom_boxplot( width=0.8, size=1) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=14),
      axis.text = element_text(size=20),
      axis.title =  element_text(size=20)
    ) +
    scale_fill_manual(values = c("white", "grey")) +
    xlab("") + ylab("Relative expression value")

    plot1
  
}


fitPlot <- function(dat, p, p1, types){
  
  temp2 <- ddply(dat, .(line, variable), summarise, avg=mean(value), sd=sd(value))
  temp2 <- temp2[temp2$variable %in% types,]
  
  
  x <- temp2$avg[temp2$line == p]; y <- temp2$avg[temp2$line == p1] 
  x_sd <- temp2$sd[temp2$line == p]; y_sd <- temp2$sd[temp2$line == p1] 
  dat <- data.frame(x, y, x_sd, y_sd, cell_type=temp2$variable[temp2$line == p])
  lim_x <- aes(xmax = x + x_sd, xmin=x - x_sd)
  lim_y <- aes(ymax = y + y_sd, ymin=y - y_sd)
  
  print(lim_y)
  
  plot1 <- ggplot(dat, aes(x, y, colour=cell_type)) + 
    geom_point(size=2) + 
    geom_errorbar(lim_y, width=0.1, size=1.2) +
    geom_errorbarh(lim_x, height=0.1, size=1.2) +
    theme_bw() + 
    theme(
      axis.text = element_text(size=15),
      axis.title =  element_text(size=15),
      legend.text = element_text(size=15)
    ) + 
    geom_abline(slope = 1, intercept = 0, lty=3) +
    xlab(p) + ylab(p1)
  
  plot1
  
}

fitPlot_1 <- function(dat, dat1, p, p1, types){

  temp1 <- ddply(dat, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
  temp2 <- ddply(dat1, .(variable), plyr::summarise, avg=mean(value), sd=sd(value))
  temp <- merge(temp1, temp2, by="variable")
  
  temp <- temp[temp$variable %in% types,]
  
  x <- temp$avg.x 
  y <- temp$avg.y
  x_sd <- temp$sd.x
  y_sd <- temp$sd.y
  dat <- data.frame(x, y, x_sd, y_sd, cell_type=temp$variable)
  lim_x <- aes(xmax = x + x_sd, xmin=x - x_sd)
  lim_y <- aes(ymax = y + y_sd, ymin=y - y_sd)
  
  plot1 <- ggplot(dat, aes(x, y, colour=cell_type)) + 
    geom_point(size=2) + 
    geom_errorbar(lim_y, width=0.1, size=1.2) +
    geom_errorbarh(lim_x, height=0.1, size=1.2) +
    theme_bw() + 
    theme(
      axis.text = element_text(size=15),
      axis.title =  element_text(size=15),
      legend.text = element_text(size=15)
    ) + 
    geom_abline(slope = 1, intercept = 0, lty=3) +
    xlab(p) + ylab(p1)
  
  plot1
  
}






## HEATMAP #############################

heatmap <- function(rep.maov){
  
  rep.maov[rep.maov == 0] <- "ns."
  rep.maov[rep.maov == -1] <- "NA"
  # rep.maov[rep.maov <= 0.01] <- "***"
  rep.maov[rep.maov <= 0.05] <- "*"
  rep.maov[rep.maov > 0.05 & rep.maov != "ns." & rep.maov != "NA"] <- "ns."
  
  dat <- as.data.frame(rep.maov)
  dat$line_1 <- rownames(dat)
  dat <- melt(dat, id.vars = c("line_1"))
  
  dat$line_2 <- dat$variable
  dat$p_value <- dat$value
  
  
  plot1 <- ggplot(dat, aes(line_1, line_2, z= p_value)) + 
    geom_tile(aes(fill = p_value)) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=15),
      axis.text = element_text(size=15),
      legend.text = element_text(size=10),
      legend.title = element_text(size=10)) +
    
    scale_fill_manual(values=c("#00BC47", "white", "grey"), 
                      name="Significance level",
                      labels=c("*", "ns.", "NA")) +
    #scale_fill_gradient(low="blue", high="white", space="Lab")   +
    xlab("") + ylab("") + 
    coord_fixed()
  
  plot1
}



heatmap_fit <- function(rep.maov, diff=F){
  
  dat <- as.data.frame(rep.maov)
  dat$line_1 <- rownames(dat)
  dat <- melt(dat, id.vars = c("line_1"))
  
  dat$line_2 <- dat$variable
  dat$p_value <- dat$value
  
  plot1 <- ggplot(dat, aes(line_1, line_2, z= p_value)) + 
    geom_tile(aes(fill = p_value)) + 
    theme_bw() + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size=12),
      axis.text = element_text(size=12),
      legend.text = element_text(size=10),
      legend.title = element_text(size=10)) +
    # scale_fill_distiller(palette="Spectral")+
    xlab("") + ylab("") + 
    coord_fixed()
  
  if(diff) plot1 <- plot1 + scale_fill_gradientn(colours=c("#EF7C09", "#FFFFFF","#51A9F9"), limits = c(-1,1), na.value = "grey")
  if(!diff) plot1 <- plot1 + scale_fill_gradientn(colours=cscale, name="r-squared values")
  
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
