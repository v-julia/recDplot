library(ggplot2)
library("ape")
library("seqinr")
library(sjstats)
library("spider")
library("gridExtra")
#library(hexbin)
library(grid)
library(gplots)
library(colorRamps)
library(stringr)


fig_path = "D:/MY_FILES/DATA/Lukashev/RNA_viruses/Lagovirus/figures/"
path_al = "D:/MY_FILES/DATA/Lukashev/RNA_viruses/Lagovirus/"

lago_name = 'Lagovirus_full_aln_0.5_less_amb_genotyped.fasta'
lago = read.dna(paste(path_al, lago_name, sep=""), format="fasta", as.character=TRUE)
lago[lago=='-'] <- NA

dna_object = lago

plot_PDCP_control(lago, 1, 3000, 3001, 7000)

#' Function plots pairwise nucleotide distance comparison plot with control.
#'
#' Each dot corresponds to a pair of nucleotide distances between
#' the same pair of genomes in two genomic regions - st1-e1 and st2-e2 (see axis).
#' Returns list with ggplot, matrices of distances between pairs of seqences calculated for
#' st1-e1 and st2-e2 regions
#' @export

plot_PDCP_control = function(dna_object, st1,e1,st2,e2){
  
# subalignments for odd and even sites
al_odd=dna_object[1:length(dna_object[,1]), seq(from = 1, to = length(dna_object[1,]), by=2)]
al_even=dna_object[1:length(dna_object[,1]), seq(from = 2, to = length(dna_object[1,]), by=2)]

# distance matrices for each region
dna_sl_dist1_odd <-dist.gene(al_odd, method = "percentage",  pairwise.deletion = TRUE)
dna_sl_dist2_even <-dist.gene(al_even, method = "percentage",  pairwise.deletion = TRUE)

# adding random noise to distances matrices' values
dist1_odd = as.vector(dna_sl_dist1_odd) + rnorm(length(dna_sl_dist1_odd),mean = 0,sd= 0.0001)
dist2_even = as.vector(dna_sl_dist2_even) + rnorm(length(dna_sl_dist2_even),mean = 0,sd= 0.0001)

dna_sl1=dna_object[1:length(dna_object[,1]), seq(from = st1, to = e1, by=1)]
dna_sl2=dna_object[1:length(dna_object[,1]), seq(from = st2, to = e2, by=1)]

# distance matrices for each region
dna_sl_dist1 <-dist.gene(dna_sl1, method = "percentage",  pairwise.deletion = TRUE)
dna_sl_dist2 <-dist.gene(dna_sl2, method = "percentage",  pairwise.deletion = TRUE)


# adding random noise to distances matrices' values

dist1= as.vector(dna_sl_dist1) + rnorm(length(dna_sl_dist1),mean = 0,sd= 0.0001)
dist2= as.vector(dna_sl_dist2) + rnorm(length(dna_sl_dist2),mean = 0,sd= 0.0001)


dist_df = data.frame(dist1_odd, dist2_even, dist1, dist2)


ggplot() +  stat_bin2d(data = dist_df, aes(dist1, dist2), binwidth = 0.003) + scale_fill_gradientn(colours=c("blue","red")) + 
  stat_bin2d(data = dist_df, aes(dist1_odd, dist2_even, alpha=0.5), binwidth = 0.003, color='black') 
#  + theme(legend.justification=c(1,0), legend.position=c(1,0))
#  + xlab(paste(toString(st1),toString(e1),sep=":"))+ylab(paste(toString(st2),toString(e2),sep=":"))


ggplot() +  geom_point(data = dist_df, aes(dist1, dist2), color="red") + 
  geom_point(data = dist_df, aes(dist1_odd, dist2_even, alpha=0.5))

}  


#g1 = ggplot(data.frame(dist1_odd,dist2_even),aes(dist1,dist2))+ stat_bin2d()+#stat_bin2d(binwidth = 0.003)+
#  scale_fill_gradientn(colours=c("blue","red"))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+
#  xlab("odd sites")+ylab("even sites")
#g2 = ggplot(data.frame(dist1,dist2),aes(dist1,dist2))+ stat_bin2d()+#stat_bin2d(binwidth = 0.003)+
#  scale_fill_gradientn(colours=c("blue","red"))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+
#  xlab("1")+ylab("2")



plot_PDCP_control = function(dna_object){
  
  # subalignments for odd and even sites
  al_odd=dna_object[1:length(dna_object[,1]), seq(from = 1, to = length(dna_object[1,]), by=2)]
  al_even=dna_object[1:length(dna_object[,1]), seq(from = 2, to = length(dna_object[1,]), by=2)]
  
  # distance matrices for each region
  dna_sl_dist1 <-dist.gene(al_odd, method = "percentage",  pairwise.deletion = TRUE)
  dna_sl_dist2 <-dist.gene(al_even, method = "percentage",  pairwise.deletion = TRUE)
  
  # adding random noise to distances matrices' values
  dist1= as.vector(dna_sl_dist1) + rnorm(length(dna_sl_dist1),mean = 0,sd= 0.0001)
  dist2= as.vector(dna_sl_dist2) + rnorm(length(dna_sl_dist2),mean = 0,sd= 0.0001)
  
  # pairwise nucleotide distance comparison plot
  dist_plot=ggplot(data.frame(dist1,dist2),aes(dist1,dist2))+ stat_bin2d()+#stat_bin2d(binwidth = 0.003)+
    scale_fill_gradientn(colours=c("blue","red"))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+
    xlab("odd sites")+ylab("even sites")
  #+  geom_smooth(method='lm',formula=y~x)
  
  return(list(dist_plot, dna_sl_dist1, dna_sl_dist2))
  
}

#' Function plots pairwise nucleotide distance comparison plot.
#'
#' Each dot corresponds to a pair of nucleotide distances between
#' the same pair of genomes in two genomic regions - st1-e1 and st2-e2 (see axis).
#' Returns list with ggplot, matrices of distances between pairs of seqences calculated for
#' st1-e1 and st2-e2 regions
#' @export
plot_dist_test = function(dna_object, st1,e1,st2,e2){
  
  # subalignments for st1-e1 and st2-e2 regions
  dna_sl1=dna_object[1:length(dna_object[,1]), seq(from = st1, to = e1, by=1)]
  dna_sl2=dna_object[1:length(dna_object[,1]), seq(from = st2, to = e2, by=1)]
  
  # distance matrices for each region
  dna_sl_dist1 <-dist.gene(dna_sl1, method = "percentage",  pairwise.deletion = TRUE)
  dna_sl_dist2 <-dist.gene(dna_sl2, method = "percentage",  pairwise.deletion = TRUE)
  #HepadnaDist1 <-dist.dna(Hepadna1, as.matrix = TRUE, model = "JC69")
  #HepadnaDist4 <-dist.dna(Hepadna4, as.matrix = TRUE, model = "JC69")
  
  # adding random noise to distances matrices' values
  
  dist1= as.vector(dna_sl_dist1) + rnorm(length(dna_sl_dist1),mean = 0,sd= 0.0001)
  dist2= as.vector(dna_sl_dist2) + rnorm(length(dna_sl_dist2),mean = 0,sd= 0.0001)
  
  
  
  # pairwise nucleotide distance comparison plot
  
  dist_plot=ggplot(data.frame(dist1,dist2),aes(dist1,dist2))+stat_bin2d(binwidth = 0.002)+
    scale_fill_gradientn(colours=c("blue","red"))+ theme(legend.justification=c(1,0), legend.position=c(1,0))+
    xlab(paste(toString(st1),toString(e1),sep=":"))+ylab(paste(toString(st2),toString(e2),sep=":"))
  #+  geom_smooth(method='lm',formula=y~x)
  
  return(list(dist_plot, dna_sl_dist1, dna_sl_dist2))
  
}
