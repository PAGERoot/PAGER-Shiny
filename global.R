
# Guillaume Lobet - University of Liege


##################################################################
###### Load libraries
##################################################################

packages <- c("ggplot2", "plyr", "gridExtra", "XML", "MASS", "reshape2", "RColorBrewer", "rsvg", "readxl", "png", "lattice", "plotly")
for(p in packages){
  if (!require(p,character.only = TRUE)){
    install.packages(p, dep=TRUE)
    if(!require(p,character.only = TRUE)) stop("Package not found")
  }
}

##################################################################
#######   Custom functions
##################################################################

normalize <- function(rs){
  for(i in 1:nrow(rs)){
    rs[i,2:ncol(rs)] <- rs[i,2:ncol(rs)] - min(rs[i,2:ncol(rs)])
    rs[i,2:ncol(rs)] <- round(rs[i,2:ncol(rs)] / max(rs[i,2:ncol(rs)]), 2) * 100
  }
  return(rs)
}

getColor <- function(val){
  co <- round(val*10)
  if(min(co) < 0) co <- co + abs(min(co))
  if(min(co) > 0) co <- co - abs(min(co))
  return(co + 1)
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}



pal <- c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')
#pal <- c("white", "red")
pal.div <- c('#b2182b','#f7f7f7','#2166ac')

