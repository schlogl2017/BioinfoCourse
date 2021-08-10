################################################################
## Installation of R and BioConductor libraries required to run the scripts of this
## course

################################################################
## R packages
install.packages(c("mda", ## For multivariate analysis
                   "class", ## Required by mda
                   "MASS", ## For multivariate analysis
                   "combinat", ## For practical on clustering
                   "maptree",  ## For practical on clustering
                   "e1071",  ## For practical on clustering
                   "ade4", ## Princomp with group-specific coloring
                   "igraph", ## Network analysis and visualization
                   "limma"
                   ))

################################################################
## Installation of the BioConductor packages required for this course.
## Last update: Feb 2011.
## source("http://bioconductor.org/biocLite.R") ## Load the R script install package
## 
BiocManager::install(c('limma', 'affy', 'arrayQuality', 'GEOquery', "hgu95av2.db", "hgu133a.db", "GO.db", 'geneplotter', 'multtest', 'golubEsets', 'ShortRead', 'siggenes'))
## Not sure I need this
BiocManager("ALL")
BiocManager('simpleaffy')

################################################################
## OBSOLETE (I don't think I still need to load the following
## libraries, to be checked).
##
## ## Packages for the 2007 version of the course
## bioconductor.packages <- c('Biobase',
##                            'affy',
##                            'simpleaffy',
##                            'marray',
##                            )
## to.install <- c('sma',
##                 'limma',
##                 bioconductor.packages)

## ## The packages I used around 2003
## install.packages("VR")
## install.packages("multiv")
## install.packages("PHYLOGR")
## install.packages("pixmap")
## install.packages("polymars")
## install.packages("mclust")
## install.packages("Matrix")
## install.packages("maptree")
## install.packages("permax")

