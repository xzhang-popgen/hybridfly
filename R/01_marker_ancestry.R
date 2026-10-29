args <- commandArgs(trailingOnly=TRUE)
data_dir <- if(length(args)>=1) args[1] else "data/marker"
out_dir <- if(length(args)>=2) args[2] else "results/marker"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

library(RColorBrewer)
library(ggplot2)
library(reshape2)

sech_maj <- read.csv(file.path(data_dir,"fly_anc-pos_sechelulans.csv"),stringsAsFactors=F)
simu_maj <- read.csv(file.path(data_dir,"fly_anc-pos_simellia.csv"),stringsAsFactors=F)

sech_maj <- sech_maj[c(1:109),]
sech_maj$order <- c(1:nrow(sech_maj))
simu_maj$order <- c(1:nrow(simu_maj))

sech_maj$chrom <- sub("_.*", "", sech_maj$Marker)
simu_maj$chrom <- sub("_.*", "", simu_maj$Marker)
chromosomes <- unique(sech_maj$chrom)

gen_cols_sech <- grep("sechelulans-N",colnames(sech_maj),value=TRUE)
gen_cols_simu <- grep("simellia-N",colnames(simu_maj),value=TRUE)
sampled_gen <- as.numeric(sub(".*-N","",gen_cols_simu))

####genome-wide ancestry
plot_these <- match(c("simellia-N2","simellia-N11","simellia-N45"),colnames(simu_maj))
symbols <- c(16,17,18)
line_type <- c(3,2,1)
line_width <- c(0.7,1,1.5)
point_size <- c(0.5,0.7,1)
ticks <- c(1,19,42,61,86,90,109)
text_pos <- c(10,31,52,74,89,100)

pdf(file.path(out_dir,"marker_genomewide_ancestry.pdf"),width=10,height=5)
par(mar=c(4.1,5.1,4.1,11),xpd=TRUE)
plot(NULL,xlim=c(0,nrow(simu_maj)),ylim=c(0,1),xlab="Genomic Positions",
     ylab="Sechellia Ancestry",xaxt="n",cex.axis=1.2,cex.lab=1.2)
axis(1,at=ticks,labels=c("","","","","","",""),tck=-0.03)
text(text_pos,-0.07,labels=c("2L","2R","3L","3R","4","X"),adj=1,xpd=TRUE,cex=1.2)
lines(x=sech_maj$order,y=sech_maj$expected,col="gray60",lwd=3,lty=2)
lines(x=simu_maj$order,y=simu_maj$expected,col="black",lwd=3,lty=2)
for(n in 1:length(plot_these)){
  points(x=sech_maj$order,y=sech_maj[,sub("simellia","sechelulans",colnames(simu_maj)[plot_these[n]])],col="red",pch=symbols[n],cex=point_size[n])
  lines(x=sech_maj$order,y=sech_maj[,sub("simellia","sechelulans",colnames(simu_maj)[plot_these[n]])],col="red",lwd=line_width[n],lty=line_type[n])
  points(x=simu_maj$order,y=simu_maj[,plot_these[n]],col="blue",pch=symbols[n],cex=point_size[n])
  lines(x=simu_maj$order,y=simu_maj[,plot_these[n]],col="blue",lwd=line_width[n],lty=line_type[n])
}
legend("topright",inset=c(-0.62,0),legend=c("simu_majority","sech_majority"),lty=1,col=c("blue","red"),bty="n",seg.len=0.5,cex=1,lwd=1.2)
legend("topright",inset=c(-0.3,0.15),legend=c("N2","N11","N45"),pch=symbols,lty=line_type,lwd=line_width,merge=T,bty="n",cex=1,seg.len=0.5)
dev.off()

####mean ancestry by chromosome
anc_time_sech_mean <- data.frame(chr=chromosomes)
anc_time_simu_mean <- data.frame(chr=chromosomes)
for(g in seq_along(gen_cols_sech)){
  anc_time_sech_mean[[as.character(sampled_gen[g])]] <- sapply(chromosomes,function(chr) mean(sech_maj[sech_maj$chrom==chr,gen_cols_sech[g]],na.rm=T))
  anc_time_simu_mean[[as.character(sampled_gen[g])]] <- sapply(chromosomes,function(chr) mean(simu_maj[simu_maj$chrom==chr,gen_cols_simu[g]],na.rm=T))
}

write.csv(anc_time_sech_mean,file.path(out_dir,"chromosome_mean_sech_majority.csv"),row.names=F)
write.csv(anc_time_simu_mean,file.path(out_dir,"chromosome_mean_simu_majority.csv"),row.names=F)

pdf(file.path(out_dir,"chromosome_ancestry_time_series.pdf"),width=8,height=5)
plot(NULL,xlim=c(1,45),ylim=c(0,1.05),xlab="Generations",ylab="Sechellia Ancestry",xaxt="n",cex.lab=1.2)
axis(side=1,at=c(1,11,27,45),paste0("N",c(1,11,27,45)),las=3,cex.axis=1)
pch_list <- c(0:4,8)
lty_list <- c(1:6)
for(n in 1:6){
  points(x=c(1,sampled_gen),y=c(0.75,as.numeric(anc_time_sech_mean[n,-1])),col="red",pch=pch_list[n],cex=0.7)
  lines(x=c(1,sampled_gen),y=c(0.75,as.numeric(anc_time_sech_mean[n,-1])),col="red",lwd=1,lty=lty_list[n])
  points(x=c(1,sampled_gen),y=c(0.25,as.numeric(anc_time_simu_mean[n,-1])),col="blue",pch=pch_list[n],cex=0.7)
  lines(x=c(1,sampled_gen),y=c(0.25,as.numeric(anc_time_simu_mean[n,-1])),col="blue",lwd=1.5,lty=lty_list[n])
}
legend("topright",inset=c(-0.55,0),legend=c("simu_majority","sech_majority"),lty=1,col=c("blue","red"),bty="n",seg.len=0.5,cex=1)
legend("topright",inset=c(-0.22,0.15),legend=chromosomes,pch=pch_list,lty=lty_list,lwd=1,merge=T,bty="n",cex=1,seg.len=0.5)
dev.off()

X_dist <- as.numeric(anc_time_simu_mean[anc_time_simu_mean$chr=="X",-1])
auto_dist <- unlist(anc_time_simu_mean[anc_time_simu_mean$chr!="X",-1])
print(wilcox.test(X_dist,auto_dist,alternative="less"))

####outlier markers
outlier <- simu_maj[simu_maj$`simellia-N45`>=0.25 & complete.cases(simu_maj[,gen_cols_simu]),]
write.csv(outlier,file.path(out_dir,"simu_majority_outlier_markers.csv"),row.names=F)

pdf(file.path(out_dir,"outlier_marker_time_series.pdf"),width=8,height=5)
plot(NULL,xlim=c(1,45),ylim=c(0,1),xlab="Generations",ylab="Sechellia Ancestry",xaxt="n",cex.lab=1.2)
axis(side=1,at=c(1,11,27,45),paste0("N",c(1,11,27,45)),las=3,cex.axis=1)
chr3R <- anc_time_simu_mean[anc_time_simu_mean$chr=="3R",-1]
chr4 <- anc_time_simu_mean[anc_time_simu_mean$chr=="4",-1]
points(x=c(1,sampled_gen),y=c(0.25,as.numeric(chr3R)),col="purple",pch=16,cex=1)
lines(x=c(1,sampled_gen),y=c(0.25,as.numeric(chr3R)),col="purple",lwd=3)
points(x=c(1,sampled_gen),y=c(0.25,as.numeric(chr4)),col="orange",pch=16,cex=1)
lines(x=c(1,sampled_gen),y=c(0.25,as.numeric(chr4)),col="orange",lwd=3)
width <- c(2,1.5,1,0.5)
for(i in seq_len(min(4,nrow(outlier)))){
  points(x=c(1,sampled_gen),y=c(0.25,as.numeric(outlier[i,gen_cols_simu])),col="purple",cex=0.8)
  lines(x=c(1,sampled_gen),y=c(0.25,as.numeric(outlier[i,gen_cols_simu])),col="purple",lwd=width[i],lty=2)
}
dev.off()

####effective selection parameter
get_s_approx <- function(p0,pn,n){
  ratio = p0*(1-pn)/(pn*(1-p0))
  s = (-log(ratio))/n
  return(s)
}

simu_freq <- simu_maj[,c("Marker","chrom",gen_cols_simu)]
for(n in 3:ncol(simu_freq)){
  simu_freq[which(simu_freq[,n]==0),n] <- 0.001
  simu_freq[which(simu_freq[,n]==1),n] <- 0.999
}

make_s <- function(g1,g2){
  c1 <- paste0("simellia-N",g1)
  c2 <- paste0("simellia-N",g2)
  get_s_approx(simu_freq[[c1]],simu_freq[[c2]],g2-g1)
}

simu_s_large <- data.frame(
  Marker=simu_freq$Marker,
  chrom=simu_freq$chrom,
  N2_11=make_s(2,11),
  N11_45=make_s(11,45),
  N2_45=make_s(2,45),
  stringsAsFactors=F
)

write.csv(simu_s_large,file.path(out_dir,"effective_selection_large_intervals.csv"),row.names=F)

simu_s_long <- melt(simu_s_large,id.vars=c("Marker","chrom"))
p <- ggplot(simu_s_long,aes(x=variable,y=Marker,fill=value))+
  geom_tile()+
  geom_text(aes(label=round(value,2)),size=2)+
  scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=0,name="m")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90,hjust=1),axis.text=element_text(color="black",size=6))+
  labs(x="Generation gaps",y="SNP Marker")
ggsave(file.path(out_dir,"effective_selection_heatmap.pdf"),p,width=6,height=10)

cat("\n5.9 Mb:\n")
print(simu_s_large[simu_s_large$Marker=="3R_05.9",])
cat("\n12.9 Mb:\n")
print(simu_s_large[simu_s_large$Marker=="3R_12.9",])
