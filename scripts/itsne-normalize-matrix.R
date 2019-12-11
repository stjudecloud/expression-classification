if(!require(optparse)) install.packages("optparse", repos = "http://cran.us.r-project.org")
if(!require(plotly)) install.packages("plotly", repos = "http://cran.us.r-project.org")
if(!require(Rtsne)) install.packages("Rtsne", repos = "http://cran.us.r-project.org")
if(!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "http://cran.us.r-project.org")
if(!require(DESeq2)) BiocManager::install("DESeq2")
if(!require(sva)) BiocManager::install("sva")

library(optparse)
library(DESeq2)
library(sva)
library(Rtsne)
library(plotly)
 
option_list = list(
  make_option(
    c("-f", "--filename"), 
    type="character", 
    default=NULL, 
    help="File containing matrix weights in CSV format.", 
    metavar="character"
  ),
  make_option(
    c("-o", "--outname"), 
    type="character", 
    default="tsne.html", 
    help="Output filename", 
    metavar="character"
  )
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

allData <- read.csv(file=opt$filename, header=T)
diagnosis <- allData["Diagnosis"][,1]
covariates <- allData["Covariates"][,1]
samples <- allData["Sample"][,1]
counts <- t(allData[!(names(allData) %in% c("Diagnosis", "Covariates", "Sample"))])
rm(allData)

colors        <- rainbow(length(unique(diagnosis)))
names(colors) <- unique(diagnosis)

# dataMatrix <- vst(as.matrix(counts), blind=T)
dataMatrix <- varianceStabilizingTransformation(as.matrix(counts), blind=T)
dataMatrix <- ComBat(dataMatrix, covariates)

topn          <- 1000
mads          <- apply(dataMatrix, 1, mad)
sortMads      <- as.numeric(sort(mads, decreasing=TRUE))
madCutoff     <- sortMads[(topn+1)]
topGenes      <- mads > madCutoff
dataMatrixTop <- dataMatrix[topGenes,]

distMat <- dist(t(dataMatrixTop))

set.seed(0)
tsne_out <- Rtsne(distMat, dims = 2, perplexity = 30, 
                  theta = 0, max_iter = 5000, check_duplicates = F )

toPlot <- data.frame(tsne_out$Y)
colnames(toPlot) <- c("xcoord","ycoord")
unknownX<-toPlot[nrow(toPlot),1]
unknownY<-toPlot[nrow(toPlot),2]
toPlot$classes   <- diagnosis

p <- plot_ly(type = "scatter" , mode = "markers" , data = toPlot[1:(nrow(toPlot)-1),], 
             x = ~xcoord, y = ~ycoord , color = ~classes , colors = colors , text = ~classes , 
             marker = list(size = 10, line = list(color = 'rgba(0,0,0)',width = 2)))%>%
  add_trace(data = toPlot[nrow(toPlot),] , x = ~xcoord, y = ~ycoord , 
            marker = list(color = 'rgb(0, 0, 0 , 0.25)', size = 20))

# Save scatterplot as interactive plotly html
htmlwidgets::saveWidget(as_widget(p), opt$outname)
