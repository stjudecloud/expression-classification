suppressMessages(library(optparse, quietly = T))
suppressMessages(library(Rtsne, quietly = T))
suppressMessages(library(plotly, quietly = T))
suppressMessages(library(DESeq2, quietly = T))
suppressMessages(library(sva, quietly = T))
suppressMessages(library("pracma", quietly = T) )
suppressMessages(library("data.table", quietly = T)) 

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
    default=20,
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
  ),
  make_option(
    c("--input-sample"),
    type="character",
    default=NULL,
    help="Sample to highlight in output graph",
    metavar="character"
  ),
  make_option(
    c("--tissue-type"),
    type="character",
    default="", 
    help="Tissue type for the input sample(s) [blood,brain,solid]",
    metavar="character"
  ),
  make_option(
    c("--save-data"),
    type="logical",
    action="store_true",
    default=F,
    help="Save 2 dimensional distance matrix from t-SNE"
  ),
  make_option(
    c("--gene-list"),
    type="character",
    default=NULL,
    help="File containing gene list for selection. See the documentation for how to format this file.",
    metavar="character"
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

# Read full data matrix
allData    <- read.csv(file=opt$filename, header=T)
# Split metadata out into separate data structures
diagnosis  <- allData["Diagnosis"][,1]
diagnosisName <- allData["DiagnosisName"][,1]
colors <- allData["Color"][,1]
projects <- allData["Projects"][,1]
covariates <- allData["Covariates"][,1]
samples    <- allData["Sample"][,1]
preservative <- allData["Preservative"][,1]
# Remove metadata from data matrix
counts     <- t(allData[!(names(allData) %in% c("Diagnosis", "Covariates", "Sample", "DiagnosisName", "Color", "Projects", "Preservative"))])

rm(allData)

# Converts to a matrix (each sample is now a column)
dataMatrix <- as.matrix(counts)
write.table(dataMatrix, file="original_matrix.txt", sep="\t",quote=FALSE,row.names=FALSE)

if (!opt$`disable-variance-stabilization`) {
  if (opt$vst) {
    cat("Performing vst\n", file = stderr())
    dataMatrix <- vst(dataMatrix, blind=T)
    write.table(dataMatrix, file="vst.txt", sep="\t",quote=FALSE,row.names=FALSE)

  } else {
    cat("Performing variance stabilization\n", file = stderr())
    dataMatrix <- varianceStabilizingTransformation(dataMatrix, blind=T)
    write.table(dataMatrix, file="variance_stabilized.txt", sep="\t",quote=FALSE,row.names=FALSE)

  }
} else {
  cat("Variance stabilization disabled!\n", file = stderr())
}

if (!opt$`disable-batch-correction`) {
  cat("Performing batch correction\n", file = stderr())
  dataMatrix <- ComBat(dataMatrix, covariates, mean.only = FALSE)
  write.table(dataMatrix, file="batch_corrected.txt", sep="\t",quote=FALSE,row.names=FALSE)

} else {
  cat("Batch correction disabled!\n", file = stderr())
}

if (is.null(opt$`gene-list`)){
  cat("Selecting gene list\n", file=stderr())
  topn          <- opt$`keep-top-n-genes`
  mads          <- apply(dataMatrix, 1, mad)
  sortMads      <- as.numeric(sort(mads, decreasing=TRUE))
  madCutoff     <- sortMads[(topn+1)]
  topGenes      <- mads > madCutoff
  dataMatrixTop <- dataMatrix[topGenes,]
} else {
  cat("Using specified gene list\n", file=stderr())
  # Read the specified gene list
  gene_list <- read.csv(file=opt$`gene-list`, header=F)
  # genes are the rows in sorted order
  genes <- row.names(counts)
  topGenes <- genes %in% gene_list$V1
  dataMatrixTop <- dataMatrix[topGenes,]
}

write.table(topGenes, file="genes.txt", sep="\t",quote=FALSE,row.names=FALSE)
write.table(dataMatrixTop, file="data_top_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)

distMat <- dist(t(dataMatrixTop))

set.seed(opt$seed)
cat("Running Rtsne\n", file = stderr())
tsne_out <- Rtsne(distMat, dims = 2, perplexity = opt$`tsne-perplexity`,
                  theta = opt$`tsne-theta`, max_iter = opt$`tsne-max-iterations`, check_duplicates = F,
                  num_threads = 0 )
cat("Saving Rtsne output\n", file = stderr())
write.table(data.frame(tsne_out$Y), file="tsne_output.txt", sep="\t",quote=FALSE,row.names=FALSE)

cat("Plotting...\n", file = stderr())
# Create plot data objects
toPlot <- data.frame(tsne_out$Y)
colnames(toPlot) <- c("t1", "t2")
unknownX<-toPlot[nrow(toPlot),1]
unknownY<-toPlot[nrow(toPlot),2]
toPlot$classes   <- diagnosis
toPlot$samples <- samples
toPlot$color <- colors
toPlot$diagnosisNames <- diagnosisName
toPlot$projects <- projects

# Setup axis
ax <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)
title = "St. Jude Cloud"

if(length(opt$`tissue-type`)){
  title = paste(title, opt$`tissue-type`)
}
title = paste(title, "RNA-Seq Expression Landscape")
plotData <- toPlot

if (opt$`save-data`) {
  cat("Saving plot data to file\n", file = stderr())
  write.table(plotData, file="tsne.txt", sep="\t",quote=FALSE,row.names=FALSE)
}

# If we have input samples, remove them from the initial plotting set and plot them separately later.
if (length(opt$`input-sample`)){
   inputs <- strsplit(opt$`input-sample`, ',')
   '%!in%' <- function(x,y)!('%in%'(x,y))
   plotData <- plotData[plotData$sample %!in% unlist(inputs) | plotData$classes %!in% unlist(inputs),]
}

p <- plot_ly()

for (category in sort(unique(plotData$classes)))
{
    subdata <- plotData[ plotData$classes %in% c(category), ]
    color <- subdata[1,5]

    p <- add_trace(p, type = "scatter" , mode = "markers" , data = subdata,
             x = ~t1, y = ~t2 , name = category,
             marker = (list( color = ~color )),
             hoverinfo = "text",
             text = ~paste("Sample: ", samples, '<br>Diagnosis Code: ', classes, '<br>Diagnosis Name: ', diagnosisNames),
             showlegend = TRUE
             )

}
p <- layout(p, title=title, xaxis=ax, yaxis=ax)

# If we have input samples, add them to the existing plot.
if (length(opt$`input-sample`)){
   inputs <- strsplit(opt$`input-sample`, ',')

   L <- toPlot[toPlot$samples %in% unlist(inputs) & toPlot$classes %in% unlist(inputs),]
   L$classes <- L$samples

   # Add input markers
   p <- add_trace(p, type = "scatter" , mode = "markers", data = L , x = ~t1, y = ~t2 , name = ~classes,
            text = ~classes, textposition = "bottom center", hoverinfo = "text", hovertext= ~classes,
            marker = list(symbol = "circle-cross",
            line = list(color="white", width=2), color="black", size=15))

   # Add input labels
   a <- list(
     x = L$t1,
     y = L$t2,
     text = L$samples,
     xref = "x",
     yref = "y",
     showarrow = TRUE,
     arrowhead = 0,
     standoff = 9,
     clicktoshow = "onoff",
     ax = 20,
     ay = -40,
     bgcolor="white",
     bordercolor="black",
     borderpad=4,
     borderwidth=2,
     arrowwidth=2
   )
   p <- layout(p, annotations=a)
}
# Save scatterplot as interactive plotly html
htmlwidgets::saveWidget(as_widget(p), opt$outname)
