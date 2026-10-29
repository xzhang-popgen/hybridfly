args <- commandArgs(trailingOnly=TRUE)
data_file <- if(length(args)>=1) args[1] else "data/readcount_by_bottle_Illumina_data"
out_dir <- if(length(args)>=2) args[2] else "results/illumina"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

library(data.table)
library(ggplot2)
library(zoo)
library(patchwork)

read_fly_counts <- function(file){
  d <- fread(file,skip="# 1-snp")
  first11 <- c(
    "snp","chr","pos",
    "simreads_sim_ref","secreads_sec_ref",
    "simreads_sec_ref","secreads_sim_ref",
    "simreads_sim_nonref","secreads_sec_nonref",
    "simreads_sec_nonref","secreads_sim_nonref"
  )
  bottle_cols <- unlist(lapply(1:61,function(b){
    c(paste0("b",b,"_nsim"),paste0("b",b,"_nsec"))
  }))
  setnames(d,c(first11,bottle_cols))
  d[,pos:=as.numeric(pos)]
  d[]
}

counts <- read_fly_counts(data_file)

#bottle groups
sim_bottles <- c(paste0("b",1:45),"b47","b60")
sec_bottles <- c("b46",paste0("b",48:58))
ref_bottles <- c("b59","b61")

#QC
counts[,total_nonref:=simreads_sim_nonref+secreads_sec_nonref+simreads_sec_nonref+secreads_sim_nonref]
counts[,parental_ok:=simreads_sim_ref>0 & secreads_sec_ref>0 & simreads_sec_ref==0 & secreads_sim_ref==0]
counts_f <- counts[parental_ok==TRUE & total_nonref>=14 & total_nonref<=44]

print(nrow(counts_f))
print(table(counts_f$chr))

make_windows <- function(d,window=300000){
  d[,bin:=floor(pos/window)*window+window/2]
  d
}

calc_pooled_windows <- function(d,window=300000){
  d <- copy(d)
  d <- make_windows(d,window)

  sim_nsim_cols <- paste0(sim_bottles,"_nsim")
  sim_nsec_cols <- paste0(sim_bottles,"_nsec")
  sec_nsim_cols <- paste0(sec_bottles,"_nsim")
  sec_nsec_cols <- paste0(sec_bottles,"_nsec")

  d[,simelia_nsim:=rowSums(.SD),.SDcols=sim_nsim_cols]
  d[,simelia_nsec:=rowSums(.SD),.SDcols=sim_nsec_cols]
  d[,sechel_nsim:=rowSums(.SD),.SDcols=sec_nsim_cols]
  d[,sechel_nsec:=rowSums(.SD),.SDcols=sec_nsec_cols]

  win <- d[,.(
    simelia_nsim=sum(simelia_nsim),
    simelia_nsec=sum(simelia_nsec),
    sechel_nsim=sum(sechel_nsim),
    sechel_nsec=sum(sechel_nsec),
    n_snps=.N
  ),by=.(chr,bin)]

  win[,simelia_p_sech:=simelia_nsec/(simelia_nsim+simelia_nsec)]
  win[,sechel_p_sech:=sechel_nsec/(sechel_nsim+sechel_nsec)]
  win[,simelia_depth:=simelia_nsim+simelia_nsec]
  win[,sechel_depth:=sechel_nsim+sechel_nsec]
  win[]
}

calc_bottle_windows <- function(d,bottles,window=300000){
  d <- copy(d)
  d[,bin:=floor(pos/window)]

  rbindlist(lapply(bottles,function(b){
    nsim <- paste0(b,"_nsim")
    nsec <- paste0(b,"_nsec")
    tmp <- d[,.(
      nsim=sum(get(nsim)),
      nsec=sum(get(nsec)),
      n_snp=.N
    ),by=.(chr,bin)]
    tmp[,bottle:=b]
    tmp[,depth:=nsim+nsec]
    tmp[,p_sech:=nsec/depth]
    tmp
  }))
}

win300 <- calc_pooled_windows(counts_f,300000)
win30 <- calc_pooled_windows(counts_f,30000)
win300[,chr:=factor(chr,levels=c("X","2L","2R","3L","3R","4"))]
win30[,chr:=factor(chr,levels=c("X","2L","2R","3L","3R","4"))]
win300[,pos_mb:=bin/1e6]
win30[,pos_mb:=bin/1e6]

fwrite(win300,file.path(out_dir,"pooled_300kb_windows.csv"))
fwrite(win30,file.path(out_dir,"pooled_30kb_windows.csv"))

####bottle-level ancestry
bw300_sim <- calc_bottle_windows(counts_f,sim_bottles,300000)
bw300_sim[,pos_mb:=bin*0.3]

chr4_b <- bw300_sim[chr=="4" & pos_mb>=0.12 & pos_mb<=1.3,
                    .(chr4=weighted.mean(p_sech,w=depth,na.rm=T)),by=bottle]
r3p1_b <- bw300_sim[chr=="3R" & pos_mb>=6.6 & pos_mb<=10.5,
                    .(r3_peak1=weighted.mean(p_sech,w=depth,na.rm=T)),by=bottle]
r3p2_b <- bw300_sim[chr=="3R" & pos_mb>=13.5 & pos_mb<=16.5,
                    .(r3_peak2=weighted.mean(p_sech,w=depth,na.rm=T)),by=bottle]
r3p3_b <- bw300_sim[chr=="3R" & pos_mb>=24.0 & pos_mb<=25.8,
                    .(r3_peak3=weighted.mean(p_sech,w=depth,na.rm=T)),by=bottle]

reg_b <- Reduce(function(x,y) merge(x,y,by="bottle"),list(chr4_b,r3p1_b,r3p2_b,r3p3_b))
fwrite(reg_b,file.path(out_dir,"bottle_region_ancestry.csv"))

print(summary(reg_b))
print(cor(reg_b[,.(chr4,r3_peak1,r3_peak2,r3_peak3)],use="complete.obs"))

####correlation plots
corr_long <- melt(reg_b,id.vars=c("bottle","chr4"),
                  measure.vars=c("r3_peak1","r3_peak2","r3_peak3"),
                  variable.name="peak",value.name="r3_ancestry")
corr_long[,peak:=factor(peak,
  levels=c("r3_peak1","r3_peak2","r3_peak3"),
  labels=c("3R peak 1","3R peak 2","3R peak 3"))]

stat_dt <- corr_long[!is.na(chr4) & !is.na(r3_ancestry),{
  fit <- lm(r3_ancestry~chr4)
  sm <- summary(fit)
  .(slope=coef(fit)[["chr4"]],pval=sm$coefficients["chr4","Pr(>|t|)"],r2=sm$r.squared)
},by=peak]
stat_dt[,label:=paste0("R2 = ",round(r2,2),"\nP = ",format(pval,digits=2,scientific=T))]
fwrite(stat_dt,file.path(out_dir,"chr4_3R_correlations.csv"))

make_corr_panel <- function(peak_name,y_label){
  dat <- corr_long[peak==peak_name]
  lab <- stat_dt[peak==peak_name]
  ggplot(dat,aes(x=chr4,y=r3_ancestry))+
    geom_smooth(method="lm",se=TRUE,linewidth=0.6)+
    geom_point(size=1.8,alpha=0.8)+
    annotate("text",x=0.05,y=0.67,hjust=0,vjust=1,label=lab$label,size=2.8)+
    coord_cartesian(xlim=c(0,1),ylim=c(0,0.7))+
    labs(x=expression("Chr4 "*italic("D. sechellia")*" ancestry"),y=y_label)+
    theme_bw(base_size=9)+theme(panel.grid.minor=element_blank())
}

p_corr1 <- make_corr_panel("3R peak 1","3R peak 1\nD. sechellia ancestry")
p_corr2 <- make_corr_panel("3R peak 2","3R peak 2\nD. sechellia ancestry")
p_corr3 <- make_corr_panel("3R peak 3","3R peak 3\nD. sechellia ancestry")
ggsave(file.path(out_dir,"chr4_vs_3R_peak1.pdf"),p_corr1,width=3.2,height=3)
ggsave(file.path(out_dir,"chr4_vs_3R_peak2.pdf"),p_corr2,width=3.2,height=3)
ggsave(file.path(out_dir,"chr4_vs_3R_peak3.pdf"),p_corr3,width=3.2,height=3)

####repeatability scan
background <- bw300_sim[
  !(chr=="4" |
    (chr=="3R" & pos_mb>=6.6 & pos_mb<=10.5) |
    (chr=="3R" & pos_mb>=13.5 & pos_mb<=16.5) |
    (chr=="3R" & pos_mb>=24.0 & pos_mb<=25.8)),
  median(p_sech,na.rm=T)
]
threshold <- background+0.10

repeat_win <- bw300_sim[,.(
  mean_p=mean(p_sech,na.rm=T),
  median_p=median(p_sech,na.rm=T),
  frac_high=mean(p_sech>threshold,na.rm=T),
  n_high=sum(p_sech>threshold,na.rm=T),
  n_bottles=sum(!is.na(p_sech)),
  mean_depth=mean(depth,na.rm=T)
),by=.(chr,bin,pos_mb)]
repeat_win[,chr:=factor(chr,levels=c("X","2L","2R","3L","3R","4"))]
repeat_win[,region:="Other"]
repeat_win[chr=="4" & pos_mb>=0.12 & pos_mb<=1.3,region:="chr4"]
repeat_win[chr=="3R" & pos_mb>=6.6 & pos_mb<=10.5,region:="3R peak1"]
repeat_win[chr=="3R" & pos_mb>=13.5 & pos_mb<=16.5,region:="3R peak2"]
repeat_win[chr=="3R" & pos_mb>=24.0 & pos_mb<=25.8,region:="3R peak3"]
fwrite(repeat_win,file.path(out_dir,"repeatability_300kb.csv"))

repeat_summary <- data.table(
  region=c("chr4","3R peak1","3R peak2","3R peak3","Background"),
  frac=c(
    mean(reg_b$chr4>threshold,na.rm=T),
    mean(reg_b$r3_peak1>threshold,na.rm=T),
    mean(reg_b$r3_peak2>threshold,na.rm=T),
    mean(reg_b$r3_peak3>threshold,na.rm=T),
    mean(bw300_sim[!(chr=="4" | (chr=="3R" & pos_mb>=6.6 & pos_mb<=10.5) | (chr=="3R" & pos_mb>=13.5 & pos_mb<=16.5) | (chr=="3R" & pos_mb>=24 & pos_mb<=25.8)),p_sech]>threshold,na.rm=T)
  )
)
fwrite(repeat_summary,file.path(out_dir,"repeatability_region_summary.csv"))

####figure 4
pA <- ggplot(repeat_win,aes(x=pos_mb,y=frac_high))+
  geom_line(linewidth=0.45)+
  facet_wrap(~chr,scales="free_x",ncol=1)+
  coord_cartesian(ylim=c(0,1))+
  geom_hline(yintercept=0.10,linetype=2)+
  labs(x="Position (Mb)",y="Fraction of bottles with elevated\nD. sechellia ancestry")+
  theme_bw(base_size=9)+theme(strip.background=element_rect(fill="grey90"),panel.grid.minor=element_blank())

reg_long <- melt(reg_b,id.vars="bottle",variable.name="region",value.name="p_sech")
reg_long[,region:=factor(region,
  levels=c("chr4","r3_peak1","r3_peak2","r3_peak3"),
  labels=c("chr4","3R peak1","3R peak2","3R peak3"))]
pB <- ggplot(reg_long,aes(region,p_sech))+
  geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.12,alpha=0.7,size=1.5)+
  coord_cartesian(ylim=c(0,1))+labs(x=NULL,y="D. sechellia ancestry")+theme_bw()

pC <- ggplot(repeat_summary,aes(reorder(region,frac),frac))+
  geom_col()+geom_text(aes(label=paste0(round(frac*100),"%")),vjust=-0.4)+
  coord_cartesian(ylim=c(0,1.05))+labs(x=NULL,y="Fraction of bottles retained")+theme_bw()

fig4 <- pA | (pB/pC)
ggsave(file.path(out_dir,"Figure4_repeatability.pdf"),fig4,width=11,height=7.5)

####30kb repeatability
bw30_sim <- calc_bottle_windows(counts_f,sim_bottles,30000)
bw30_sim[,pos_mb:=bin*0.03]
background30 <- bw30_sim[
  !(chr=="4" |
    (chr=="3R" & pos_mb>=6.6 & pos_mb<=10.5) |
    (chr=="3R" & pos_mb>=13.5 & pos_mb<=16.5) |
    (chr=="3R" & pos_mb>=24.0 & pos_mb<=25.8)),
  median(p_sech,na.rm=T)
]
threshold30 <- background30+0.10
repeat30 <- bw30_sim[,.(
  mean_p=mean(p_sech,na.rm=T),
  median_p=median(p_sech,na.rm=T),
  frac_high=mean(p_sech>threshold30,na.rm=T),
  n_high=sum(p_sech>threshold30,na.rm=T),
  n_bottles=sum(!is.na(p_sech)),
  mean_depth=mean(depth,na.rm=T)
),by=.(chr,bin,pos_mb)]
fwrite(repeat30,file.path(out_dir,"repeatability_30kb.csv"))

####figure 3
p3A <- ggplot(win300,aes(x=pos_mb,y=simelia_p_sech))+
  geom_line(linewidth=0.45)+facet_wrap(~chr,scales="free_x",ncol=1)+
  coord_cartesian(ylim=c(0,1))+
  labs(x="Position (Mb)",y=expression(italic("D. sechellia")~"ancestry"))+
  theme_bw(base_size=9)+theme(strip.background=element_rect(fill="grey90"),panel.grid.minor=element_blank())

chr4_30 <- win30[chr=="4"]
p3B <- ggplot(chr4_30,aes(x=pos_mb,y=simelia_p_sech))+
  geom_line(linewidth=0.6)+coord_cartesian(ylim=c(0,1))+
  labs(x="Chromosome 4 position (Mb)",y=expression(italic("D. sechellia")~"ancestry"))+
  theme_bw(base_size=9)+theme(panel.grid.minor=element_blank())

r3_30 <- win30[chr=="3R"]
p3C <- ggplot(r3_30,aes(x=pos_mb,y=simelia_p_sech))+
  geom_line(linewidth=0.5)+coord_cartesian(ylim=c(0,0.65),xlim=c(0,28))+
  annotate("rect",xmin=6.6,xmax=10.5,ymin=-Inf,ymax=Inf,alpha=0.10)+
  annotate("rect",xmin=13.5,xmax=16.5,ymin=-Inf,ymax=Inf,alpha=0.10)+
  annotate("rect",xmin=24.0,xmax=25.8,ymin=-Inf,ymax=Inf,alpha=0.10)+
  labs(x="Chromosome 3R position (Mb)",y=expression(italic("D. sechellia")~"ancestry"))+
  theme_bw(base_size=9)+theme(panel.grid.minor=element_blank())

peak_summary <- data.table(
  region=c("Genome background","3R peak 3","3R peak 1","3R peak 2","chr4"),
  mean_ancestry=c(background,
    mean(reg_b$r3_peak3,na.rm=T),
    mean(reg_b$r3_peak1,na.rm=T),
    mean(reg_b$r3_peak2,na.rm=T),
    mean(reg_b$chr4,na.rm=T))
)
peak_summary[,region:=factor(region,levels=c("Genome background","3R peak 3","3R peak 1","3R peak 2","chr4"))]
p3D <- ggplot(peak_summary,aes(x=region,y=mean_ancestry))+
  geom_col(width=0.65)+coord_cartesian(ylim=c(0,1))+
  labs(x=NULL,y=expression("Mean "*italic("D. sechellia")~"ancestry"))+
  theme_bw(base_size=9)+theme(panel.grid.minor=element_blank(),axis.text.x=element_text(angle=30,hjust=1))

ggsave(file.path(out_dir,"Fig3A_genomewide_300kb.pdf"),p3A,width=6,height=6)
ggsave(file.path(out_dir,"Fig3B_chr4_zoom_30kb.pdf"),p3B,width=4.5,height=3)
ggsave(file.path(out_dir,"Fig3C_3R_zoom_30kb.pdf"),p3C,width=6,height=3)
ggsave(file.path(out_dir,"Fig3D_region_mean_ancestry.pdf"),p3D,width=4.5,height=3.2)

cat("background ancestry:",background,"\n")
cat("retention threshold:",threshold,"\n")
print(repeat_summary)
