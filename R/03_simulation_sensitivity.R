args <- commandArgs(trailingOnly=TRUE)
sim_dir <- if(length(args)>=1) args[1] else "results/simulations"
out_dir <- if(length(args)>=2) args[2] else "results/simulation_summary"
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

library(data.table)
library(ggplot2)

make_pseudo <- function(dat,n_exp=10000,n_bottles=46){
  reps <- unique(dat$rep)
  rbindlist(lapply(seq_len(n_exp),function(i){
    chosen <- sample(reps,n_bottles,replace=T)
    sample_table <- data.table(sample_id=seq_along(chosen),rep=chosen)
    x <- merge(sample_table,dat,by="rep",allow.cartesian=T)
    x[,.(f59=mean(f_5.9),f129=mean(f_12.9)),by=generation][,experiment:=i]
  }))
}

####neutral drift
neutral <- fread(file.path(sim_dir,"neutral_1000.csv"))
neutral_long <- melt(neutral,
  id.vars=c("rep","generation"),
  measure.vars=c("f_5.9","f_12.9"),
  variable.name="locus",value.name="freq")
neutral_long[,locus:=sub("f_","",locus)]
neutral_traj <- neutral_long[,.(
  mean=mean(freq),
  q025=quantile(freq,0.025),
  q975=quantile(freq,0.975)
),by=.(generation,locus)]
fwrite(neutral_traj,file.path(out_dir,"neutral_trajectory.csv"))

p_neutral <- ggplot(neutral_traj,aes(generation,mean))+
  geom_ribbon(aes(ymin=q025,ymax=q975),alpha=0.2)+
  geom_line()+facet_wrap(~locus)+theme_classic()+
  labs(x="Generation",y="D. sechellia ancestry frequency")
ggsave(file.path(out_dir,"neutral_trajectory.pdf"),p_neutral,width=7,height=4)

set.seed(12345)
pseudo_neutral <- make_pseudo(neutral,10000,46)
pseudo_neutral_stats <- pseudo_neutral[,{
  p2_59 <- f59[generation==2]
  p11_59 <- f59[generation==11]
  p45_59 <- f59[generation==45]
  p2_129 <- f129[generation==2]
  p11_129 <- f129[generation==11]
  p45_129 <- f129[generation==45]
  .(
    early59=p11_59-p2_59,late59=p45_59-p11_59,
    early129=p11_129-p2_129,late129=p45_129-p11_129,
    empirical59=(p11_59-p2_59)<=-0.23 & (p45_59-p11_59)>=0.45,
    empirical129=(p11_129-p2_129)<=-0.20 & (p45_129-p11_129)>=0.23
  )
},by=experiment]

neutral_summary <- pseudo_neutral_stats[,.(
  mean_early59=mean(early59),mean_late59=mean(late59),
  p_empirical59=mean(empirical59),
  mean_early129=mean(early129),mean_late129=mean(late129),
  p_empirical129=mean(empirical129)
)]
fwrite(neutral_summary,file.path(out_dir,"neutral_46bottle_summary.csv"))

####recombination
rec_files <- c(
  "0"="AE20_rec_0_500.csv",
  "0.25"="AE20_rec_0_25_500.csv",
  "0.5"="AE20_rec_0_5_500.csv",
  "1"="AE20_rec_1_500.csv",
  "2"="AE20_rec_2_500.csv"
)
all_rec <- rbindlist(lapply(names(rec_files),function(x){
  d <- fread(file.path(sim_dir,rec_files[[x]]))
  d[,rec_mult:=as.numeric(x)]
  d
}),fill=T)

rec_summary <- all_rec[,.(
  p0_59=mean(f_5.9[generation==0]),
  p11_59=mean(f_5.9[generation==11]),
  p45_59=mean(f_5.9[generation==45]),
  p0_129=mean(f_12.9[generation==0]),
  p11_129=mean(f_12.9[generation==11]),
  p45_129=mean(f_12.9[generation==45])
),by=rec_mult]
rec_summary[,`:=`(
  early59=p11_59-p0_59,late59=p45_59-p11_59,
  early129=p11_129-p0_129,late129=p45_129-p11_129
)]
fwrite(rec_summary,file.path(out_dir,"recombination_summary.csv"))

bottle_rec <- all_rec[,{
  p11_59 <- f_5.9[generation==11]
  p45_59 <- f_5.9[generation==45]
  p11_129 <- f_12.9[generation==11]
  p45_129 <- f_12.9[generation==45]
  .(
    late59=p45_59-p11_59,rebound05_59=p45_59>=p11_59+0.05,
    late129=p45_129-p11_129,rebound05_129=p45_129>=p11_129+0.05
  )
},by=.(rec_mult,rep)]

recomb_prob <- bottle_rec[,.(
  rebound05_59=mean(rebound05_59),
  rebound05_129=mean(rebound05_129)
),by=rec_mult]
fwrite(recomb_prob,file.path(out_dir,"recombination_bottle_rebound.csv"))

p_recomb <- ggplot(melt(recomb_prob,id.vars="rec_mult",variable.name="locus",value.name="prob"),
                   aes(rec_mult,prob,linetype=locus))+
  geom_line()+geom_point()+theme_classic()+
  labs(x="Recombination rate multiplier",y="P(G45 >= G11 + 0.05)",linetype="Locus")
ggsave(file.path(out_dir,"recombination_rebound.pdf"),p_recomb,width=6,height=4)

set.seed(12345)
pseudo_rec <- rbindlist(lapply(sort(unique(all_rec$rec_mult)),function(r){
  p <- make_pseudo(all_rec[rec_mult==r],10000,46)
  p[,rec_mult:=r]
  p
}))
pseudo_rec_stats <- pseudo_rec[,{
  p11_59 <- f59[generation==11]
  p45_59 <- f59[generation==45]
  p11_129 <- f129[generation==11]
  p45_129 <- f129[generation==45]
  .(
    late59=p45_59-p11_59,reverse59=p45_59>p11_59,
    late129=p45_129-p11_129,reverse129=p45_129>p11_129
  )
},by=.(rec_mult,experiment)]
pseudo_rec_summary <- pseudo_rec_stats[,.(
  mean_late59=mean(late59),reverse59=mean(reverse59),
  mean_late129=mean(late129),reverse129=mean(reverse129)
),by=rec_mult]
fwrite(pseudo_rec_summary,file.path(out_dir,"recombination_46bottle_summary.csv"))

####population size
N10 <- fread(file.path(sim_dir,"AE20_N10_500.csv")); N10[,Ne:=10]
N20 <- fread(file.path(sim_dir,"AE20_N20_500.csv")); N20[,Ne:=20]
N30 <- fread(file.path(sim_dir,"AE20_N30_500.csv")); N30[,Ne:=30]
N50 <- fread(file.path(sim_dir,"AE20_rec_1_500.csv")); N50[,Ne:=50]
allN <- rbindlist(list(N10,N20,N30,N50),fill=T)

N_summary <- allN[,.(
  p0_59=mean(f_5.9[generation==0]),p11_59=mean(f_5.9[generation==11]),p45_59=mean(f_5.9[generation==45]),
  p0_129=mean(f_12.9[generation==0]),p11_129=mean(f_12.9[generation==11]),p45_129=mean(f_12.9[generation==45])
),by=Ne]
N_summary[,`:=`(
  early59=p11_59-p0_59,late59=p45_59-p11_59,
  early129=p11_129-p0_129,late129=p45_129-p11_129
)]
fwrite(N_summary,file.path(out_dir,"population_size_summary.csv"))

bottleN <- allN[,{
  p11_59 <- f_5.9[generation==11]; p45_59 <- f_5.9[generation==45]
  p11_129 <- f_12.9[generation==11]; p45_129 <- f_12.9[generation==45]
  .(
    survive45_59=p45_59>0,rebound05_59=p45_59>=p11_59+0.05,
    survive45_129=p45_129>0,rebound05_129=p45_129>=p11_129+0.05
  )
},by=.(Ne,rep)]
N_prob <- bottleN[,.(
  survive45_59=mean(survive45_59),rebound05_59=mean(rebound05_59),
  survive45_129=mean(survive45_129),rebound05_129=mean(rebound05_129)
),by=Ne]
fwrite(N_prob,file.path(out_dir,"population_size_bottle_summary.csv"))

N_plot <- rbind(
  N_prob[,.(Ne,locus="5.9",prob=rebound05_59)],
  N_prob[,.(Ne,locus="12.9",prob=rebound05_129)]
)
p_N <- ggplot(N_plot,aes(Ne,prob,linetype=locus))+
  geom_line()+geom_point()+theme_classic()+
  labs(x="Population size",y="P(G45 >= G11 + 0.05)",linetype="Locus")
ggsave(file.path(out_dir,"population_size_rebound.pdf"),p_N,width=6,height=4)

set.seed(12345)
pseudoN <- rbindlist(lapply(c(10,20,30,50),function(n){
  p <- make_pseudo(allN[Ne==n],10000,46)
  p[,Ne:=n]
  p
}))
pseudoN_stats <- pseudoN[,{
  p11_59 <- f59[generation==11]; p45_59 <- f59[generation==45]
  p11_129 <- f129[generation==11]; p45_129 <- f129[generation==45]
  .(
    late59=p45_59-p11_59,reverse59=p45_59>p11_59,rebound05_59=p45_59>=p11_59+0.05,
    late129=p45_129-p11_129,reverse129=p45_129>p11_129,rebound05_129=p45_129>=p11_129+0.05
  )
},by=.(Ne,experiment)]
pseudoN_summary <- pseudoN_stats[,.(
  mean_late59=mean(late59),reverse59=mean(reverse59),rebound05_59=mean(rebound05_59),
  mean_late129=mean(late129),reverse129=mean(reverse129),rebound05_129=mean(rebound05_129)
),by=Ne]
fwrite(pseudoN_summary,file.path(out_dir,"population_size_46bottle_summary.csv"))

####historical sampling
sampled <- rbind(
  N50[sample_n>0 & obs_5.9>=0,.(rep,generation,sample_n,locus="5.9",true_freq=f_5.9,obs_freq=obs_5.9)],
  N50[sample_n>0 & obs_12.9>=0,.(rep,generation,sample_n,locus="12.9",true_freq=f_12.9,obs_freq=obs_12.9)]
)
sampled[,error:=obs_freq-true_freq]
sampling_error <- sampled[,.(
  bias=mean(error),MAE=mean(abs(error)),RMSE=sqrt(mean(error^2))
),by=.(sample_n,locus)]
fwrite(sampling_error,file.path(out_dir,"sampling_per_bottle_error.csv"))

make_sampling_pseudo <- function(dat,n_exp=10000,n_bottles=46){
  reps <- unique(dat$rep)
  rbindlist(lapply(seq_len(n_exp),function(i){
    chosen <- sample(reps,n_bottles,replace=T)
    sample_table <- data.table(sample_id=seq_along(chosen),rep=chosen)
    x <- merge(sample_table,dat,by="rep",allow.cartesian=T)
    x[,.(true_mean=mean(true_freq),observed_mean=mean(obs_freq)),by=.(generation,sample_n,locus)][,experiment:=i]
  }))
}

set.seed(12345)
pseudo_sampling <- make_sampling_pseudo(sampled,10000,46)
pseudo_sampling[,error:=observed_mean-true_mean]
pseudo_sampling_summary <- pseudo_sampling[,.(
  bias=mean(error),MAE=mean(abs(error)),RMSE=sqrt(mean(error^2)),
  q025=quantile(error,0.025),q975=quantile(error,0.975)
),by=.(sample_n,locus)]
fwrite(pseudo_sampling_summary,file.path(out_dir,"sampling_46bottle_error.csv"))

phase_accuracy <- pseudo_sampling[generation %in% c(2,10,44),{
  true2 <- true_mean[generation==2]; true10 <- true_mean[generation==10]; true44 <- true_mean[generation==44]
  obs2 <- observed_mean[generation==2]; obs10 <- observed_mean[generation==10]; obs44 <- observed_mean[generation==44]
  true_early <- true10-true2; obs_early <- obs10-obs2
  true_late <- true44-true10; obs_late <- obs44-obs10
  .(
    true_early=true_early,obs_early=obs_early,true_late=true_late,obs_late=obs_late,
    early_ok=sign(true_early)==sign(obs_early),late_ok=sign(true_late)==sign(obs_late),
    both_ok=sign(true_early)==sign(obs_early) & sign(true_late)==sign(obs_late)
  )
},by=.(experiment,locus)]
phase_summary <- phase_accuracy[,.(
  mean_true_early=mean(true_early),mean_observed_early=mean(obs_early),early_direction_accuracy=mean(early_ok),
  mean_true_late=mean(true_late),mean_observed_late=mean(obs_late),late_direction_accuracy=mean(late_ok),
  two_phase_accuracy=mean(both_ok),early_MAE=mean(abs(obs_early-true_early)),late_MAE=mean(abs(obs_late-true_late))
),by=locus]
fwrite(phase_summary,file.path(out_dir,"sampling_phase_accuracy.csv"))

####linkage architecture
arch_dir <- file.path(sim_dir,"architecture")
if(dir.exists(arch_dir)){
  arch_files <- list.files(arch_dir,pattern="\\.csv$",full.names=T)
  if(length(arch_files)>0){
    arch <- rbindlist(lapply(arch_files,function(f){
      d <- fread(f)
      d[,file:=basename(f)]
      d
    }),fill=T)
    arch_summary <- arch[,.(
      p0_59=mean(f_5.9[generation==0]),p11_59=mean(f_5.9[generation==11]),p45_59=mean(f_5.9[generation==45]),
      p0_129=mean(f_12.9[generation==0]),p11_129=mean(f_12.9[generation==11]),p45_129=mean(f_12.9[generation==45])
    ),by=.(scenario,file)]
    arch_summary[,`:=`(
      early59=p11_59-p0_59,late59=p45_59-p11_59,
      early129=p11_129-p0_129,late129=p45_129-p11_129
    )]
    fwrite(arch_summary,file.path(out_dir,"linkage_architecture_summary.csv"))
  }
}

cat("\nneutral\n")
print(neutral_summary)
cat("\nrecombination\n")
print(rec_summary)
print(pseudo_rec_summary)
cat("\npopulation size\n")
print(N_summary)
print(pseudoN_summary)
cat("\nsampling\n")
print(pseudo_sampling_summary)
print(phase_summary)
