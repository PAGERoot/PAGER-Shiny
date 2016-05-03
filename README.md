# PAGE-Root | Shiny

## Overview

PAGE-Root is a tool aiming at comparing and analysing gene expression patterns in different cell layers. The gene data can be either coming from gene reporter fluoresce data, such as the ones coming from CellSet, or gene expression data coming from transcriptomic experiments. 

PAGE-Root was designed to work with root related data, but could be easily adapter to any type of spatially register data. 

## Before running PAGE-Root

Some dependencies are needed before running PAGE-Root. These are mainly needed for the display of the different root cell layers. 


### rsvg

Informations about the `rsvg` package can be found here: [https://github.com/jeroenooms/rsvg#installation]()

Binary packages for **OS-X** or **Windows** can be installed directly from CRAN:

```r
install.packages("rsvg")
```

Installation from source on Linux or OSX requires [`librsvg2`](https://developer.gnome.org/rsvg/). On **Debian** or **Ubuntu** install [librsvg2-dev](https://packages.debian.org/testing/librsvg2-dev):

```
sudo apt-get install -y librsvg2-dev
```

On **Fedora**, **CentOS or RHEL** we need [librsvg-devel](https://apps.fedoraproject.org/packages/librsvg2-devel):

```
sudo yum install librsvg-dev
```

On **OS-X** use [rsvg](https://github.com/Homebrew/homebrew-core/blob/master/Formula/librsvg.rb) from Homebrew:

```
brew install librsvg
```


### XML


For Linux users, some dependencies needs to be installed to use the `XML` package in R:

```
sudo apt-get install libcurl4-openssl-dev
sudo apt-get install libxml2-dev
```




## Running PAGE-Root

PAGE-Root has been developed as a Shiny app, to ease its maintenance and deployment. 

Running the app is really simple an does not require any coding skills:

1. Install [R](https://www.r-project.org/) (if you do not have the latest version of R, it will help to update it)
2. (*Optional*) Install [RStudio](https://www.rstudio.com/) 
3. Open `R` and, in the console, type in the following code:

``` 
install.packages("shiny")
library("shiny")
shiny::runGitHub("PAGERoot/PAGER-Shiny", "PAGERoot")
```

4. PAGE-Root will open, either in a new window if you are using `RStudio`, or in your favorite web browser.


> The first time PAGE-Root is launched, it might take some time, since it has to install multiple libraries on your system. Do not worry and take a coffee :)


## Using PAGE-Root

### Without your own dataset

- check `Use example data`
- click `Load and analyse data`
- wait for the loading to be over
- click `Plot data`
