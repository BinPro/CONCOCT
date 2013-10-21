#!/usr/bin/Rscript

#load libraries
library(ggplot2)
library(getopt)
library(grid)

spec = matrix(c('verbose','v',0,"logical",'help','h',0,"logical",'bicfile','b',1,"character",'ofile','o',1,"character",'cstart','s',2,"integer",'cend','e',2,"integer"),byrow=TRUE,ncol=4)

opt=getopt(spec)

# if help was asked for print a friendly message 
# and exit with a non-zero error code 
if( !is.null(opt$help)) {
	cat(getopt(spec, usage=TRUE)); 
	q(status=1);
}

bicFile <- opt$bicfile

BIC <- read.csv(bicFile,header=FALSE)

BIC_r <-range(BIC$V1)

xStart <- BIC_r[[1]]
xEnd <- BIC_r[[2]]

if( !is.null(opt$cstart)) {
	xStart <- opt$cstart
}

if( !is.null(opt$cend)) {
	xEnd <- opt$cend
}

BICS <- subset(BIC,BIC$V1 >= xStart & BIC$V1 <= xEnd)

BICS_y <- range(BICS$V2)

pdf(opt$ofile)

qplot(data=BICS,V1,V2) +xlim(xStart, xEnd) + ylim(BICS_y[[1]],BICS_y[[2]]) + geom_line() + ylab("BIC") + xlab("Number of components K")

dev.off()