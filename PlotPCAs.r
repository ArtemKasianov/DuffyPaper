#plot PCA in R
library(calibrate)

setwd("/directory/with/input/files") #set working directory 

evec<- read.table("finalPCA_table.05_21_25.txt", header=TRUE, sep="\t") # txt(tab delimited table) file with PCA data and colors for populations in last column
evec<- read.table("finalPCA_table.05_21_25.grey.txt", header=TRUE, sep="\t") # txt file with PCA data and colors for populations in last column(version highlighting specific populations)
eval<- read.table ("plink.eigenval") #file with eigenvalues 
namer <- "PCA_Figure_1" #Figure name
namer <- "PCA_Figure_1B" #Second figure name
nrpc<-10 

## calculate loadings for principal components
nreval<-nrow(eval)
totalev <-sum(eval)

aa <- array(NA,dim=c(nreval,1))
for (i in 1:nreval) 
  aa[i,1]<-format(round(((eval[i,1]/totalev)*100),3), nsmall = 3)}
##



pdf(file =paste(namer, "_PCA.pdf", sep=""))

# plot PCA based on data in 4 and 5 columns

grid(nx = NULL, ny = NULL,
     lty = 0,      # Grid line type
     col = "grey93", # Grid line color
     lwd = 1)
plot(-evec[,4], -evec[,5], col = evec$col, pch = evec$Pch, xlab=paste("PC1: ", aa[1,1], "%", sep=""), ylab=paste("PC2: ", aa[2,1], "%", sep=""),lwd=1,panel.first = grid(lty = 1,col = "grey93"))

axis(1, col = "gray", col.axis = "gray")  # Bottom (x-axis)
axis(2, col = "gray", col.axis = "gray")  # Left (y-axis)

box(col = "gray")
dev.off()



#Legend in separate file 

pdf(file='PCA_Legend.pdf', sep="")

par(xpd=TRUE)
plot(1,1,axes=F,xlab=" ",ylab=" ",col="white")
legend(
  "top",
  c("ǁAni-Khwe","non-ǁAni-Khwe","Kwangali","Mbukushu","Ovimbundu","Nyaneka","TaaEast","TaaNorth","TaaWest","ǂKhomani","JuǀʼhoanNorth","JuǀʼhoanSouth","Nama","Shua","Haiǁom","Somali","Datog","Masai"),
  pch=c(5,4,11,10,21,23,25,24,23,22,24,25,11,2,0,4,8,6),
  col=c("green","green","red","red","red","red","blue","blue","blue","blue","darkblue","darkblue","green","green","green","orange","orange","orange"),
  cex=0.5)
dev.off()