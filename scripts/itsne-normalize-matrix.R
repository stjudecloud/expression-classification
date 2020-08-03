suppressMessages(library(optparse, quietly = T))
suppressMessages(library(Rtsne, quietly = T))
suppressMessages(library(plotly, quietly = T))
suppressMessages(library(DESeq2, quietly = T))
suppressMessages(library(sva, quietly = T))
suppressMessages(library("RColorBrewer", quietly = T))
suppressMessages(library("viridis", quietly = T) )
suppressMessages(library("stringr", quietly = T) )
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
    cat("Performing vst\n", file = stderr())
    dataMatrix <- vst(dataMatrix, blind=T)
  } else {
    cat("Performing variance stabilization\n", file = stderr())
    dataMatrix <- varianceStabilizingTransformation(dataMatrix, blind=T)
  }
} else {
  cat("Variance stabilization disabled!\n", file = stderr())
}

if (!opt$`disable-batch-correction`) {
  cat("Performing batch correction\n", file = stderr())
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
cat("Running Rtsne\n", file = stderr())
tsne_out <- Rtsne(distMat, dims = 2, perplexity = opt$`tsne-perplexity`,
                  theta = opt$`tsne-theta`, max_iter = opt$`tsne-max-iterations`, check_duplicates = F )

# Setup colors for diagnosis codes
popcolor_All =
  c(
    "ALAL"="#d3d3d3",
    "ALCL"="#f9779d",
    "AMKL"="#008cff",
    "AMLM7"="#008cff",
    "AML"="#00c0ff",
    "AMLCBFA2T3BCL7A"="#00c0ff",#"#0000ff",
    "AMLCBFA2T3GLIS2"="#00c0ff",#"#0000ff",
    "AMLCBFBMYH11"="#00c0ff",#"#0000ff",
    "AMLDEKNUP214"="#00c0ff",
    "AMLKMT2AELL"="#00c0ff",#"#5610d1",
    "AMLKMT2AMLLT10"="#00c0ff",#"#5610d1",
    "AMLKMT2AMLLT3"="#00c0ff",#"#5610d1",
    "AMLKMT2AMLLT4"="#00c0ff",#"#5610d1",
    "AMLKMT2AMLLT6"="#00c0ff",#"#5610d1",
    "AULKMT2AMLLT3"="#00c0ff",#"#5610d1",
    "TLLKMT2AMLLT1"="red",
    "AMLM1KMT2AMLLT3"="#00c0ff",#"#5610d1",
    "AMLM5KMT2AMLLT3"="#00c0ff",#"#5610d1",
    "AMLM0KMT2AELL"="#00c0ff",#"#5610d1",
    "AMLM5KMT2AMLLT10"="#00c0ff",#"#5610d1",
    "AMLM5KMT2AMLLT1"="#00c0ff",#"#5610d1",
    "AMLM5KMT2AMLLT4"="#00c0ff",#"#5610d1",
    "AMLRUNX1CBFA2T3"="#00c0ff",#"#0000ff",
    "AMLRUNX1EVX1"="#00c0ff",#"#0000ff",
    "AMLRUNX1RUNX1T1"="#00c0ff",#"#0000ff",
    "APLPMLRARA"="#00c0ff",#"#ffa500",
    "BL"="#d3d3d3",
    "B-ALL"="grey",
    "BLLBCRABL1"="#ff00ff",
    "Ph"="#ff00ff",
    "Ph-like"="#9759d5",
    "BLLBCRABL1JAK2"="#9759d5",
    "BLLBCRABL1L"="#9759d5",
    "BLLBCRABL1LABL1"="#9759d5",
    "BLLBCRABL1LCRLF2"="#9759d5",
    "BLLBCRABL1LEPOR"="#9759d5",
    "BLLDUX4"="#696969",
    "DUX4"="#696969",
    "BLLETV6RUNX1"="#ffd700",
    "ETV6"="#ffd700",
    "BLLETV6RUNX1L"="#c6af34",
    "BLLHYPER"="#3E9F32",
    "Hyperdiploid"="#3E9F32",
    "BLLHYPO"="#483d8b",
    "Hypodiploid"="#483d8b",
    "BLLIAMP21"="#0000ff",
    "iAMP21"="#0000ff",
    "BLLIGHCEBPD"="#d3d3d3",
    "BLLKMT2A"="#7cfc00",
    "KMT2A"="#7cfc00",

    
    "BLLKMT2AUSP2"="#7cfc00",
    "BLLKMT2AMLLT10"="#7cfc00",
    "BLLKMT2AAFF1"="#7cfc00",
    "BLLKMT2AMLLT3"="#7cfc00",
    "BLLKMT2AMLLT1"="#7cfc00",
    "BLLKMT2APAK2"="#7cfc00",
    "BLLKMT2ARELA"="#7cfc00",
    
    "BLLKMT2Ainf"="#11c333",
    "BLLMEF2D"="#66C2A6",
    "BLLMYC"="#d3d3d3",
    "BLLNOS"="#d3d3d3",
    "BLLNUTM1"="#8b0000",
    "NUTM1"="#8b0000",
    "BLLPAX5ALT"="#e88c38",
    "BLLPAX5P80R"="#ffa500",
    "PAX5alt"="#e88c38",
    "PAX5 P80R"="#ffa500",
    "BLLTCF3HLF"="#daa520",
    "BLLTCF3PBX1"="#c8a2c8",
    "TCF3-PBX1"="#c8a2c8",
    "BLLZNF384"="#A8DD00",
    "ZNF384"="#A8DD00",
    "BLLZNF384L"="#8fb90a",
    "CML"="#90ee90",
    "CUP"="#d3d3d3",#"black",
    "Unknown Hematological"="#d3d3d3",
    "Unknown Solid Tumor"="#d3d3d3",
    "ST"="#d3d3d3",
    "DLBCLNOS"="#d3d3d3",
    "LCH"="#d3d3d3",
    "MDS"="#d3d3d3",
    "MS"="#d3d3d3",
    "Rare Hematological"="#d3d3d3",
    "TLL"="red",
    "T-ALL"="red",
    "TALL"="red",
    "ACC"="#66C2A6",


    "ASPS"="#d3d3d3",

    "CCRCC"="#d3d3d3",
    "DES"="#8b0000",
    "Rare Solid Tumor"="#8b0000",
    
    "DSRCT"="#daa520",
    
    #solid germ cell tumors
    "DYSNOS"="#ffd700",#"yellow",
    "GCT"="#ffd700",#"yellow",
    "MGCT"="#ffd700",#"yellow",
    "ODYS"="#ffd700",#"yellow",
    "OMGCT"="#ffd700",#"yellow",
    "YSTNOS"="#ffd700",#"yellow",
 
    "ARMS"="#00aeff", 
    "ERMS"="#0000ff",
    "BERMS"="#a3d7ff",#"#0006c2",
    "SCRMS"="#d3d3d3",#"#a3d7ff",
    "RMS"="#d3d3d3",#"#00c0ff",
    "RHB"="#d3d3d3",#"#00c0ff",   

    "ES"="#d277f3",
    "GIST"="#d3d3d3",
    
    "HCC"="#e88c38",#"#ffa500",
    "LM"="#e88c38",
    "Liver Malignancy"="#e88c38",
    "LIHB"="#e88c38",#"#e76836",
    "UESL"="#e88c38",
    
    "MEL"="#9531ed",


    "FIB"="#d3d3d3",
    "MUCC"="#d3d3d3",
    "IFS"="#d3d3d3",
    "LGFMS"="#d3d3d3",
    "NFIB"="#d3d3d3",
    "SIPT"="#d3d3d3",
    "MFH"="#d3d3d3",
    "AFH"="#d3d3d3",
    "RCSNOS"="#d3d3d3",
    "SCCNOS"="#d3d3d3",
    "SCSNOS"="#d3d3d3",
    "SETTLE"="#d3d3d3",
    "DFSP"="#d3d3d3",
    "CCSK"="#d3d3d3",
    "CHDM"="#d3d3d3",
    "SYNS"="#d3d3d3",
    "PANET"="#d3d3d3",
    "CHOS"="#d3d3d3",
    
    "NBL"="#f9779d",
    "OS"="#ff00ff",
    

    "Rare Brain Tumor"="#eb1414",
    "PRCC"="#eb1414",#"#ff130f", #renal, papillary
    "RCC"="#eb1414", #renal
    "MRT"="#c01111",
    "MRTL"="#c01111",
    
    "RBL"="#e76836",#"#ffd700", #retinoblastoma

    "THFO"="#0ea05c",
    "THPA"="#11c598",

    "WT"="#29a20b",
    "WTB"="#7cfc00",
    "WTL"="#29a20b",
    "WTR"="#29a20b",
    "WLM"="#29a20b",
   
    #BRAIN Tumors
    
    "ACPG"="red",
    "APE"="#ffccff",
    "EPM"="#ffccff",
    "EPMT"="#ffccff",
    "EPMTPF"="#ff00ff",#"#ff00ff",
    "APEPF"="#ff00ff",
    "EPMTST"="#c042ff",
    "APEST"="#c042ff",
    "APESPT"="#ef62b6",
    "MPE"="#ffccff",#"black",
    "MEPMST"="#ffccff",#"black",
    
    "ATRT"="#f9779d",
    "BGCT"="#ff7b29",
    "BMGCT"="#ff7b29",
    "BYST"="#ff7b29",
    "EBMT"="#ff7b29",
    "ETMR"="#ff7b29",
    "CPC"="#ffd700",
  
    "GB"="#d3d3d3",
    "GNG"="#11c333",

    "DMG"="#053fff",#"#0006c2",
    "HGG"="#053fff",
    "HGGNBS"="#053fff",
    "HGGNOS"="#053fff",#"#00c0ff",

    "LGGNOS"="#00c0ff",#"yellow",
    "LGG"="#00c0ff",#"yellow",

    
    "HGNET"="#8fb90a",
    
    "MBL"="#7cfc00",
    "MBLG3"="#7cfc00",#"#2fd090",
    "MBLG4"="#7cfc00",#"#2fd090",
    "MBLSHH"="#29a20b",
    "MBLWNT"="#287415",
    "MBT"="#d3d3d3",
    "MNG"="#8b0000",

    "MPEFV"="#d3d3d3",
    "MPEPF"="#d3d3d3",
    "MPNST"="#e88c38",
    "PBL"="#696969",
    "XPA"="#d3d3d3"
  )
color_lib <- as.data.frame(popcolor_All)
setDT(color_lib, keep.rownames = TRUE)[]
colnames(color_lib) <- c('classes', 'color')
cat("Plotting...\n", file = stderr())
# Create plot data objects
toPlot <- data.frame(tsne_out$Y)
colnames(toPlot) <- c("t1", "t2")
unknownX<-toPlot[nrow(toPlot),1]
unknownY<-toPlot[nrow(toPlot),2]
toPlot$classes   <- diagnosis
toPlot$samples <- samples
toPlot <- merge(toPlot, color_lib, by="classes")
#toPlot$color <- popcolor_All[toPlot$classes]

# Setup axis
ax <- list(
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)
title = "tSNE St. Jude Cloud"

if(length(opt$`tissue-type`)){
  if(opt$`tissue-type` == "solid"){
    title = paste(title, "Solid Tumors")
  }
  else if (opt$`tissue-type` == "blood"){
    title = paste(title, "Blood Cancers")
  }
  else if (opt$`tissue-type` == "brain"){
    title = paste(title, "Brain Cancers")
  }
}

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

# Plot the reference samples
p <- plot_ly(type = "scatter" , mode = "markers" , data = plotData[1:(nrow(plotData)),],
             x = ~t1, y = ~t2 , color = ~classes , colors = popcolor_All , 
             hoverinfo = "text",
             text = ~paste("Sample: ", samples, '<br>Diagnosis: ', classes) #~samples 
             )%>%
     layout(title=title, xaxis=ax, yaxis=ax)

# If we have input samples, add them to the existing plot.
if (length(opt$`input-sample`)){
   inputs <- strsplit(opt$`input-sample`, ',')

   L <- toPlot[toPlot$samples %in% unlist(inputs) & toPlot$classes %in% unlist(inputs),]
   L$classes <- L$samples
   highlight      <- rainbow(length(unique(L$samples)))
   names(highlight) <- unique(L$samples)
   p <- add_trace(p, data = L , x = ~t1, y = ~t2 , color = ~classes, colors = highlight,
            marker = list(size = 10,color = 'rgb(0, 0, 0 , 0.25)', symbol = "cross"))
   a <- list(
     x = L$t1,
     y = L$t2,
     text = L$samples,
     xref = "x",
     yref = "y",
     showarrow = TRUE,
     arrowhead = 7,
     ax = 20,
     ay = -40
   )
   p <- layout(p, annotations=a)
}
# Save scatterplot as interactive plotly html
htmlwidgets::saveWidget(as_widget(p), opt$outname)
