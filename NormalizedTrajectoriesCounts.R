library(MASS)
library(ggplot2)
library(reshape2)
library(stringr)
numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

column = c('c1','c2','c4','c5','c7','c8','c10','c11')
plate = c('p3','p4')
# column = c('c1')
# plate = c('p3')

for (a1 in column){
  for (a2 in plate){
    input = paste(c("Trajectories/",a2,"_",a1,"_normTrajectory.txt"),collapse = '')
    dat=read.table(input,sep="\t",header=T)
    output = paste(c("Trajectories/",a2,"_",a1,"_normTrajectory.pdf"),collapse = '')
    pdf(output,width=10, height=6)
    par(mar=c(4,5,1,1))
    #print(dat)
    mdat <- melt(dat, id=c("BC"))
    #print(mdat)
    mdat$variable <- as.numeric(numextract(mdat$variable))
    print(ggplot(mdat, aes(log(value))) + geom_density(aes(fill=factor(variable)), alpha=0.8) + labs(title="Density plots per generation", subtitle="BC counts by library", x="log(counts)",y="barcodes",fill="generations"))
    print(ggplot(mdat, aes(value)) + geom_density(aes(fill=factor(variable)), alpha=0.8) + labs(title="Density plots per generation", subtitle="BC counts by library", x="counts",y="barcodes",fill="generations"))
    print(ggplot(mdat,aes(x=variable,y=log10(value),group=BC)) + geom_line(alpha=0.3) + xlab("time (generations)")+ylab("log10(normalized count)") + theme(axis.text=element_text(size=20),axis.title=element_text(size=22)) + theme(legend.text=element_text(size=4))+ theme(panel.background = element_rect(fill ="transparent", colour = "black"), panel.grid.minor = element_blank(), panel.grid.major = element_blank()))
    dev.off()
  }
}
