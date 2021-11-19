library(MASS)
library(ggplot2)
mat = c('1432','E4')

for (m in mat){
  input = paste(c("Barcodes/",m,"_AggregateCounts.txt"),collapse = '')
  dat=read.table(input,sep="\t",header=T)
  output = paste(c("Barcodes/",m,"_AggregateCounts.pdf"),collapse = '')
  pdf(output,width=10, height=6)
  par(mar=c(4,5,1,1))
  truehist(dat$libraries, main="number of libraries",cex.lab=1.75,cex.axis=1.75)
  truehist(dat$counts, main="counts",cex.lab=1.75,cex.axis=1.75)
  truehist(dat$counts, main="counts", ylim = c(0,0.0001), cex.lab=1.75,cex.axis=1.75)
  #truehist(dat$counts, main="counts", xlim = c(0,1000), ylim = c(0,0.00005), cex.lab=1.75,cex.axis=1.75)
  # truehist(dat$V2, main="Counts", xlim = c(0,5000), ylim = c(0,0.00001), cex.lab=1.75,cex.axis=1.75)
  # truehist(dat$V2, main="Counts", xlim = c(0,1000), ylim = c(0,0.00001), cex.lab=1.75,cex.axis=1.75)
  plot(dat$libraries,dat$counts)
  plot(dat$libraries,log(dat$counts))
  print(ggplot(dat, aes(log(counts))) + geom_density(aes(fill=factor(libraries)), alpha=0.8) + labs(title="Density plot", 
                                                                                              subtitle="BC counts by library",
                                                                                              x="log(counts)",
                                                                                              fill="# of libraries"))
  dev.off()
}
