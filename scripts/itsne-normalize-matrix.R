library(optparse, quietly = T)
library(Rtsne, quietly = T)
library(plotly, quietly = T)
library(DESeq2, quietly = T)
library(sva, quietly = T)
 
option_list = list(
  make_option(
   c("-f", "--filename"), 
    type="character", 
    default=NULL, 
    help="File containing matrix weights in CSV format. See the documentation for how to format this file.", 
    metavar="character"
  ),
  make_option(
    c("-o", "--outname"), 
    type="character", 
    default=NULL,
    help="Filename of the output HTML file.", 
    metavar="character"
  ),
  make_option(
    c("-v", "--vst"),
    type="logical",
    action="store_true",
    default=F,
    help="Use the vst() function instead of the varianceStabilizingTransform(). This runs much more quickly but is less robust."
  ),
  make_option(
    c("-n", "--keep-top-n-genes"),
    type="integer",
    default=1000,
    help="Keep the top `n` differentially expressed genes for the t-SNE transformation."
  ),
  make_option(
    c("--seed"),
    type="integer",
    default=0,
    help="Random seed."
  ),
  make_option(
    c("--disable-variance-stabilization"),
    action="store_true",
    default=F,
    help="Disable variance stabilization (not recommended, mostly for testing purposes)."
  ),
  make_option(
    c("--disable-batch-correction"),
    action="store_true",
    default=F,
    help="Disable batch correction (not recommended, mostly for testing purposes)."
  ),
  make_option(
    c("--tsne-perplexity"),
    type="integer",
    default=30,
    help="t-SNE perplexity parameter."
  ),
  make_option(
    c("--tsne-theta"),
    type="integer",
    default=0,
    help="t-SNE theta parameter."
  ),
  make_option(
    c("--tsne-max-iterations"),
    type="integer",
    default=5000,
    help="t-SNE maximum iterations."
  )
); 
 
opt_parser = OptionParser(prog = "itsne-normalize-matrix.R", 
                          description = "Normalize a RNA-seq counts matrix and produce a t-SNE plot.",
                          option_list = option_list);
opt = parse_args(opt_parser);

if(is.null(opt) || is.null(opt$outname) || is.null(opt$filename)) {
  print_help(opt_parser)
  q(status=1)
}

allData    <- read.csv(file=opt$filename, header=T)
diagnosis  <- allData["Diagnosis"][,1]
covariates <- allData["Covariates"][,1]
samples    <- allData["Sample"][,1]
counts     <- t(allData[!(names(allData) %in% c("Diagnosis", "Covariates", "Sample"))])
rm(allData)

colors        <- rainbow(length(unique(diagnosis)))
names(colors) <- unique(diagnosis)

dataMatrix <- as.matrix(counts)
if (!opt$`disable-variance-stabilization`) {
  if (opt$vst) {
    dataMatrix <- vst(dataMatrix, blind=T)
  } else {
    dataMatrix <- varianceStabilizingTransformation(dataMatrix, blind=T)
  }
} else {
  cat("Variance stabilization disabled!\n", file = stderr())
}

if (!opt$`disable-batch-correction`) {
  dataMatrix <- ComBat(dataMatrix, covariates)
} else {
  cat("Batch correction disabled!\n", file = stderr())
}

topn          <- opt$`keep-top-n-genes`
mads          <- apply(dataMatrix, 1, mad)
sortMads      <- as.numeric(sort(mads, decreasing=TRUE))
madCutoff     <- sortMads[(topn+1)]
topGenes      <- mads > madCutoff
dataMatrixTop <- dataMatrix[topGenes,]

distMat <- dist(t(dataMatrixTop))

set.seed(opt$seed)
tsne_out <- Rtsne(distMat, dims = 2, perplexity = opt$`tsne-perplexity`, 
                  theta = opt$`tsne-theta`, max_iter = opt$`tsne-max-iterations`, check_duplicates = F )

toPlot <- data.frame(tsne_out$Y)
colnames(toPlot) <- c("xcoord", "ycoord")
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
