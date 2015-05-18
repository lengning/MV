
options <- commandArgs(trailingOnly = TRUE)
print(options)
File=options[1] # file name
MeanLow0=options[2] # lower limit of high mean genes
MeanHigh0=options[3] # upper limit of high mean genes
Norm=options[4] # whether perform normalization; if "T" is specified, median-by-ratio normalization will be performed.
#csv=options[5] # whether output csv

if(length(options)<2)MeanLow0 <- 100
if(length(options)<3)MeanHigh0 <- "NULL"
if(length(options)<4)Norm <- "F"
#if(length(options)<5)Norm <- "T"


if(MeanLow0=="NULL")MeanLow <- NULL
if(MeanLow0!="NULL")MeanLow=as.numeric(MeanLow0)

if(MeanHigh0=="NULL")MeanHigh <- NULL
if(MeanHigh0!="NULL")MeanHigh=as.numeric(MeanHigh0)




# csv or txt
tmp=strsplit(File, split="\\.")[[1]]
FileType=tmp[length(tmp)]

if(FileType=="csv"){
	cat("\n Read in csv file \n")
	prefix=strsplit(File,split="\\.csv")[[1]][1]
	In=read.csv(File,stringsAsFactors=F,row.names=1)
}
if(FileType!="csv"){
	cat("\n Read in tab delimited file \n")
	prefix=strsplit(File,split=paste0("\\.",FileType))[[1]][1]
	In=read.table(File,stringsAsFactors=F,row.names=1, sep="\t",header=T)
}



Matraw=data.matrix(In)

cat(paste("\n",
					"high expressers are defines as the ones with mean greater than",MeanLow0,
				 "and less than", MeanHigh0, "\n"))	

Mat=Matraw
print(str(Mat))

library(Oscope)

Sizes=rep(NA,ncol(Mat))
if(Norm=="T"){
cat("\n ==== Performing normalization ==== \n")
library(EBSeq)
Sizes=MedianNorm(Mat)
if(is.na(Sizes[1]))cat("\n Warning: all genes have 0(s), normalization is not performed \n")
}

#
library(Oscope)
pdf(paste0(prefix,"_MV.pdf"))
if(Norm=="T"&& (!is.na(Sizes[1])))MV=CalcMV(Data=Mat, Sizes=Sizes, NormData=FALSE, MeanCutLow=MeanLow
																				,MeanCutHigh=MeanHigh)
if(is.na(Sizes[1])) MV=CalcMV(Data=Mat, NormData=TRUE, MeanCutLow=MeanLow
					                             ,MeanCutHigh=MeanHigh)
dev.off()

if(length(MV$GeneToUse)<1)stop("no gene passed!")

M1=cbind(MV$Mean, MV$Median, MV$Var)
colnames(M1)=c("Mean","Median","Variance")
M2=matrix(MV$GeneToUse,ncol=1)
colnames(M2)="gene_name"
#write.csv(M1,file=paste0(prefix,"_MeanMedVar.csv"))
#write.csv(M2,file=paste0(prefix,"_HighMHighV.csv"))
write.table(M1,file=paste0(prefix,"_MeanMedVar.txt"))
write.table(M2,file=paste0(prefix,"_HighMHighV.txt"))
write.table(Mat[MV$GeneToUse,],paste0(prefix,"_expression_HighMHighV.txt"))



