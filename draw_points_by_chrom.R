library(SNPchip)
library(oligoClasses)
chrlen <- getSequenceLengths("hg19")[c(paste("chr",1:22,sep=''),"chrX","chrY")]
chroms <- c(paste("chr",1:22,sep=""),"chrX","chrY")
rt <- read.table("NIPT.gc",header = F,sep="\t")
chromsindex = 1:length(chroms)

for (chrindex in chromsindex){
  chrpoints =rt[rt$V1 == chroms[chrindex],]
  pic <- paste0(chroms[chrindex],".png")
  if (length(chrpoints$V1) >2){
    print(pic)
    png(pic, width = 960,height=360)
    layout(matrix(c(1,1,2),3,1,byrow = TRUE))
    
    par(bty="n")
    plot(chrpoints$V3,chrpoints$V7*100, xlim=c(0,chrlen[chrindex]),pch=19,col='blue',main=chroms[chrindex],cex.main=2, cex=1.5,
         xaxt="n", xlab="", ylab="GC%", yaxp=c(20, 60, 5), ylim=c(20, 60), cex.axis=2)
    #segments(chrpoints$V2, chrpoints$V4 / 10, chrpoints$V3, chrpoints$V4 / 10, col="red", cex=1.5)
    z=loess(chrpoints$V7*100~chrpoints$V2)
    lines(chrpoints$V2,z$fit,cex=7.5,col="red");
    in_n <- 1000000
    in_n2 <- 1000000
    inters <- seq(0,chrlen[chrindex],by=in_n)
    inters_c <- paste0(as.character(inters /in_n2),"M")
    axis(1,inters,inters_c,cex.axis=2)
    axis(2,c(1,2,3,4,5),as.character(c(1,2,3,4,5)),cex.axis=2)
    plot(chrpoints$V3,chrpoints$V4,xlim=c(0,chrlen[chrindex]),xlab='',ylab='',axes=F, col="white")
    usr1 <- par()$usr
    usr1[3]=0
    usr1[4]=5
    par(usr=usr1,omd=c(0,1,0.052,1))
    plotIdiogram(chroms[chrindex] , build="hg19",cex=1.3, new=F, xlim=c(0, chrlen[chrindex]), ylim=c(0, 6), axes=F)
    dev.off()
  }
}