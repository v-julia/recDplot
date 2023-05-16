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
library(dplyr)

fig_path = "D:/MY_FILES/DATA/Lukashev/RNA_viruses/Lagovirus/figures/"
path_al = "D:/MY_FILES/DATA/Lukashev/RNA_viruses/Lagovirus/"

path_al = "D:/Virology/RNA_viruses/Lagovirus/"

lago_name = 'Lagovirus_full_aln_0.5_less_amb_genotyped.fasta'
lago = read.dna(paste(path_al, lago_name, sep=""), format="fasta", as.character=TRUE)
lago[lago=='-'] <- NA

dna_object = lago

dist_df_lago_list = plot_PDCP_control(lago, 1, 3000, 3001, 7000)

df1 = data.frame(dist_df_lago_list$dist1, dist_df_lago_list$dist2, factor="real")
colnames(df1) = c('dist1', 'dist2')
df2 =  data.frame(dist_df_lago_list$dist1_odd, dist_df_lago_list$dist2_even, factor="control")
colnames(df2) = c('dist1', 'dist2')
d = rbind(df1, df2)
colnames(d) = c('dist1', 'dist2', 'cond')

ggplot(data = d, aes(dist1, dist2)) +  geom_point(aes(color=cond)) +  scale_color_manual(values = c("dodgerblue4", "dimgrey"))+
  #geom_point(data = d, aes(dist1, dist2), alpha=0.5, color="dimgrey")+
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  xlab(paste(toString(1),toString(3000),sep=":"))+ylab(paste(toString(3001),toString(7000),sep=":"))+
  guides(color = guide_legend(title = ""))

ggplot() +  stat_bin2d(data = dist_df_lago, aes(dist1, dist2), binwidth = 0.003) + scale_fill_gradientn(colours=c("blue","red")) +
  stat_bin2d(data = dist_df_lago, aes(dist1_odd, dist2_even, alpha=0.5), binwidth = 0.003, color='black') +
   theme(legend.justification=c(1,0), legend.position=c(1,0)) +
   xlab(paste(toString(1),toString(3000),sep=":"))+ylab(paste(toString(3001),toString(7000),sep=":"))

ggplot() +  geom_point(data = d, aes(dist1, dist2), color="dodgerblue4") +
  geom_point(data = d, aes(dist1, dist2), alpha=0.5, color="dimgrey")+
  theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  xlab(paste(toString(1),toString(3000),sep=":"))+ylab(paste(toString(3001),toString(7000),sep=":"))


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
  stat_bin2d(data = dist_df, aes(dist1_odd, dist2_even, alpha=0.5), binwidth = 0.003, color='black') +
   theme(legend.justification=c(1,0), legend.position=c(1,0)) +
  xlab(paste(toString(st1),toString(e1),sep=":"))+ylab(paste(toString(st2),toString(e2),sep=":"))


ggplot() +  geom_point(data = dist_df, aes(dist1, dist2), color="red") +
  geom_point(data = dist_df, aes(dist1_odd, dist2_even, alpha=0.5))

return(dist_df)
}

#' plots heatmap with RMSE in pairwise distance comparison plot for each pair of genomic regions
#' dna_object -  list of DNA sequences (class DNAbin)
#' step
#' window - length of genomic regions to compare
#' method - method of calculation distances ("pdist", "JC", "Kimura", "TN")
#' modification - pairwise deletion of positions with gaps or not

#'returns matrix with rmse values for each pair f=of genomic regions
#' @export
plot_rmse_control = function(dna_object, step, window, method, modification=NA){

  length_aln = length(dna_object[1,]) #length of alignment
  num_seq = length(dna_object[,1]) # number of sequences in alignment

  starts = seq(from=0, to=length_aln-window, by = step) # start positions of genomic regions
  starts[1]=1
  ends = seq(from=window, to = length_aln, by = step) # end positions of genomic regions
  if (length_aln%%step>step){ends=c(ends,length_aln)}

  df_intervals = cbind(starts,ends) #intervals

  #names = apply(df_intervals, 1, function(x){paste(toString(x[1]),toString(x[2]),sep="_")})

  #dataframe to store RMSE values of each comparison
  rmse_df = data.frame(matrix(ncol=length(starts), nrow = length(starts)))
  colnames(rmse_df)=starts
  rownames(rmse_df)=starts


  cat("Total number of intervals", nrow(df_intervals))
  #list of distance matrices for each pair of genomic regions
  dist_matrices = list()
  for (i in 1:nrow(df_intervals)){
    cat("\r", "Calculating distances in interval", i)
    slice = dna_object[1:num_seq, seq(from = df_intervals[i,"starts"], to = df_intervals[i,"ends"], by=1)]

    #dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")
    if (method == "pdist"){
      if (modification=="pairwise"){
        dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = TRUE)}
      else {
        dist_matrices[[i]] = dist.gene(slice, method = "percentage",  pairwise.deletion = FALSE)}
    }

    else {
      if (method == "JC"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "JC69")}
      if (method == "Kimura"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      if (method == "TN"){dist_matrices[[i]] = dist.dna(slice,  as.matrix = TRUE,  model = "K80")}
      #else{print("Unknown method")}
    }

  }
  n = nrow(df_intervals)
  print(length(dist_matrices))
  for (i in 1:n){
    cat("\r", "Calculating rmse in row", i)
    for (j in i:n){
      #for (j in 1:(n)){
      #print(paste(toString(i), toString(j), sep=","))
      #fits pairwise distance comparison plots linear model, calculates rmse
      rmse_i_j = (rmse(lm(dist_matrices[[j]]~dist_matrices[[i]])) + rmse(lm(dist_matrices[[i]]~dist_matrices[[j]]))) /2.0
      #rmse_i_j = rmse(lm(dist_matrices[[j]]~dist_matrices[[i]]))
      rmse_df[i,j] = rmse_i_j
      rmse_df[j,i] = rmse_i_j
    }
  }
  #print(rmse_df)
  #colnames(rmse_df)
  return(rmse_df)

}

