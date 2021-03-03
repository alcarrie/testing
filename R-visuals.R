# Installing the required packages and dependencies, uncomment if not using AWS
#install.packages("BiocManager")
#install.packages("forcats")
#install.packages("stringr")
#BiocManager::install("GEOquery")
#BiocManager::install("limma")
#BiocManager::install("pheatmap")
#BiocManager::install("plsgenomics")

# Install GEOquery
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("GEOquery")

# Loading the required packages and dependencies
library("BiocManager")
library("forcats")
library("stringr")
library("limma")

#library(dpylr)
library(GEOquery)
library(pheatmap)
library(plsgenomics)


get_file_name <- function(filename) {

  gse <- getGEO(filename=filename, destdir= "")    # change my_id to the required dataset
  #gse <- getGEO(filename='GDS859.soft.gz', destdir=".")
  # check how many platforms used
  length(gse)

  # if more than one dataset is present, can analyse the other dataset by changing the number inside the [[...]]
  # e.g. gds <- gds[[2]]

  # Turning GDS object into an expression set object (using base 2 logarithms) and examining it
  eset <- GDS2eSet(gse, do.log2=TRUE)


  sampleNames(eset)
  featureData(eset)

  # Check the normalisation and the scales used
  ## exprs get the expression levels as a data frame and get the distribution
  summary(exprs(eset))



  # Inspect the clinical variables
  sampleInfo <- pData(eset)
  # Visualisation of Data and Exploratory Analyses
  # Sample clustering and Principal Components Analysis

  # Heatmap
  corMatrix <- cor(exprs(eset),use="c")

  # Print the rownames of the sample information and check it matches the correlation matrix
  rownames(sampleInfo)
  colnames(corMatrix)

  # Ensure rownames match the columns
  rownames(sampleInfo) <- colnames(corMatrix)

  jpeg(filename= 'static/visuals/pheatmap.jpeg', units = "px", width = 4000, height = 2451, res = 300)
  pheatmap(corMatrix,annotation_col=sampleInfo)
  dev.off()
  
  ###heatmap for transcription factors activity###
  
  #set the variable names (rownames) of X to the gene names, get rid of 2 extra columns
  X <- Table(gse)
  geneNames <- as.character(X$IDENTIFIER)
  
  X <- exprs(eset)
  rownames(X) <- geneNames
  
  #to make the results easier to interpret, we make an average of the genes that have multiple duplicates
  #we use avereps() to make the average of the probes
  #installed the limma library to use avereps()
  Z<- avereps(X)
  
  #get rid of the rows containing empty values
  Z<- Z[complete.cases(Z), ]  
  
  # keep only the 100 genes with the highest SD
  #Y2 will include only 100 genes with higher SD across all samples
  install.packages("matrixStats")     
  library(matrixStats)
  
  # Add column with standard deviation for each gene
  Y2<- Z
  Y2<- transform(Y2, SD=apply(Y2,1, sd, na.rm = TRUE))
  
  # Sort rows by SD size
  Y2<- Y2[with(Y2, order(-SD)),]
  
  # take top 100 genes
  Y2<- head(Y2, 100)
  
  # remove the SD column
  Y2<- Y2[1:(length(Y2)-1)]
  
  #converting table to a matrix
  best100 <- as.matrix(Y2)
  
  Y3 <- t(best100)
  
  #Produce a heatmap only using the 100 genes with the highest SD  
  # transpose 
  library("gplots")
  
  Y3 <- data.matrix(Y3)
  par(mar = c(0,0,0,0))
  jpeg(filename= 'static/visuals/heatmap.jpeg', units = "px", width = 4000, height = 2451, res = 300)
  heatmap.2(best100, scale="column", col=heat.colors(10), Colv=TRUE)
  dev.off()
  
  
  
  
  Y <- t(Z)
 
  Mean <- colMeans(Y)
  Zscore<- scale(Mean)
  pvalue = 2*pnorm(abs(Zscore), lower.tail = F)
  
  

  
  
  
  
  

  
  
  
  
}

error_fix <- function() {
  dev.off()
}

# for RA results
relative_activity <- function() {

  connec <- read.csv(file = 'connec_data.csv', header = FALSE)
  connec2 <- read.csv(file = 'connec_data.csv')
  connec <- data.matrix(connec, rownames.force = NA)
  connec <- as.matrix(connec)
  ge <- read.csv(file = 'ge_data.csv', header = FALSE)
  ge2 <- read.csv(file = 'ge_data.csv')
  ge <- data.matrix(ge, rownames.force = NA)
  ge <- as.matrix(ge)

  new <- TFA.estimate(CONNECdata = connec, GEdata = ge,ncomp=3,nruncv=0)
  TFAc <- new$TFA
  colnames(TFAc) <- colnames(ge2)
  rownames(TFAc) <- colnames(connec2)
  TFAc[is.nan(TFAc)] = 0
  TFAc <- TFAc[,-1]
  TFAc <- TFAc[-1,]
  #TFAc <- TFAc[ order(rowMeans(TFAc), decreasing = T), ] #doesn't work with datatables
  # Add a new column for the average relative activity
  TFAc <- cbind(TFAc, rowMeans(TFAc))
  # Rename column as average
  colnames(TFAc)[ncol(TFAc)] <- "AVERAGE"
  write.csv(TFAc, "relative.csv")

}

clear_env <- function() {
    rm(list = ls()) # clear the R environment

}
