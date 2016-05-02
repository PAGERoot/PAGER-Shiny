# PAGE-Root | Shiny

## Overview

PAGE-Root is a tool aiming at comparing and analysing gene expression patterns in different cell layers. The gene data can be either coming from gene reporter fluoresce data, such as the ones coming from CellSet, or gene expression data coming from transcriptomic experiments. 

PAGE-Root was designed to work with root related data, but could be easily adapter to any type of spatially register data. 


## Running PAGE-Root

PAGE-Root has been developed as a Shiny app, to ease its maintenance and deployment. 

Running the app is really simple an does not require any coding skills:

1. Install [R](https://www.r-project.org/)
2. (*Optional*) Install [RStudio](https://www.rstudio.com/) 
3. Open `R` and, in the console, type in the following code:
 
		install.packages("shiny")
		library("shiny")
		shiny::runGitHub("PAGERoot/PAGER-Shiny", "PAGERoot")
4. PAGE-Root will open, either in a new window if you are using `RStudio`, or in your favorite web browser.


> The first time PAGE-Root is launched, it might take some time, since it has to install multiple libraries on your system. 


## Using PAGE-Root

PAGE-Root requires two types of input from the users: 

- Fluorescence data from reporter lines
- Gene expression data