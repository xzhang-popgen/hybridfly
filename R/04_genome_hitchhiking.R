args <- commandArgs(trailingOnly=TRUE)
sim_dir <- if(length(args)>=1) args[1] else "results/simulations/genome_hitchhiking"
out_dir <- if(length(args)>=2) args[2] else "results/genome_hitchhiking_summary"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

library(data.table)
library(ggplot2)

files <- list.files(sim_dir,pattern="\\.csv$",full.names=T)
if(length(files)==0) stop("no simulation files found")
d <- rbindlist(lapply(files,fread),fill=T)

summary_dt <- d[,.(
  p0=mean(genome_mean_ancestry[generation==0]),
  p11=mean(genome_mean_ancestry[generation==11]),
  p45=mean(genome_mean_ancestry[generation==45]),
  sd45=sd(genome_mean_ancestry[generation==45])
),by=.(N,n_selected,s_del)]
summary_dt[,`:=`(early=p11-p0,late=p45-p11)]
fwrite(summary_dt,file.path(out_dir,"genome_hitchhiking_summary.csv"))

traj <- d[,.(
  mean_ancestry=mean(genome_mean_ancestry),
  q025=quantile(genome_mean_ancestry,0.025),
  q975=quantile(genome_mean_ancestry,0.975)
),by=.(N,n_selected,s_del,generation)]

p <- ggplot(traj,aes(generation,mean_ancestry,group=factor(n_selected),linetype=factor(n_selected)))+
  geom_line()+facet_grid(N~s_del)+theme_classic()+
  labs(x="Generation",y="Mean minority ancestry",linetype="Selected loci")
ggsave(file.path(out_dir,"genome_hitchhiking_trajectories.pdf"),p,width=8,height=7)

print(summary_dt)
