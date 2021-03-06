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


##################################################################
###### Load libraries
##################################################################

# packages <- c("ggplot2", "plyr", "gridExtra", "MASS",
#               "reshape2", "RColorBrewer", "readxl",
#               "lattice", "plotly", "gridExtra", "xml2", "data.table")
# for(p in packages){
#   if (!require(p,character.only = TRUE)){
#     install.packages(p, dep=TRUE)
#     if(!require(p,character.only = TRUE)) stop("Package not found")
#   }
# }
library("ggplot2")
library("plyr")
library("gridExtra")
library("MASS")
library("reshape2")
library("RColorBrewer")
library("readxl")
library("lattice")
library("plotly")
library("gridExtra")
library("xml2")
library("data.table")
library("Hmisc")


##################################################################
###### Options
##################################################################
cscale <- c("#ffffcc", "#a1dab4","#41b6c4", "#2c7fb8", "#253494")

options(shiny.trace=F)

DFLT_action_enable_scrolling <- TRUE
DFLT_scrolling_y_limit <- 200

##################################################################
#######   Custom functions
##################################################################

range01 <- function(x){
  (x-min(x))/(max(x)-min(x))
}

# Function to read the RSML file
read_rsml <- function(path){
  
  # path <- "~/Downloads/reporter.rsml"
  
  # t1 <- proc.time()
  x <- read_xml(path)
  scene <- xml_children(x)[2]
  mydata <- list()#data.frame("line"=factor(), 
                       # "root"=numeric(), 
                       # "cell_type"=factor(), 
                       # "value"=numeric())
                       # 
  # all <- xml_find_all(scene, ".//root")
  # names <- xml_attr(xml_parent(all), "label")
  # test <- NULL
  
  inc <- 0
  for(i in c(1:length(xml_children(scene)))){
    plant <- xml_child(scene, i)
    # message(paste0("importing ",xml_attr(plant, "label")))
    for(j in c(1:length(xml_children(plant)))){
      root <- xml_child(plant, j)
      annots <- xml_find_all(root, ".//annotation")
      for(k in c(1:length(annots))){
        n <- length(xml_find_all(annots[k], ".//value"))
        # mydata <- rbind(mydata,
        mydata[[length(mydata)+1]] <- data.frame(line = rep(xml_attr(plant, "label"), n),
                                           root = as.numeric(rep(xml_attr(root, "label"), n)),
                                           cell_type = rep(xml_attr(annots[k], "name"), n),
                                           value = as.numeric(xml_text(xml_find_all(annots[k], ".//value"))))#)
        # tot <- tot+n
        # inc <- inc+1
        # test <- rbind(test, data.table(inc=inc, n=n, time=(proc.time()-t1)[[3]]))
      }
    }
  }
  mydata <- ldply(mydata, data.frame)
  # print(proc.time() - t1)
  return(mydata)
}

cell_types <- c("columella","cortex","endodermis","epidermis","lateralrootcap","QC","stele" )

