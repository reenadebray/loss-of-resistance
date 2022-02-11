# Analysis of bacterial genomes evolved in the absence of phage
# Reena Debray
# Feb 1, 2022

### Functions
## This function uses the output from the program breseq to construct a matrix of genes that were mutated in each population during experimental evolution
### It takes data in the following form: "breseq", a data frame with an entry for each mutation in each population at each time and its corresponding allele frequency
### "gene_list", a list of all genes observed in the population or subpopulation of interest
### "sample_list", a list of all samples in the population or subpopulation of interest 
### AF_min, a minimum allele frequency for consideration. I set AF_min=0 (no minimum), though note that breseq only returns polymorphisms with a population frequency of 0.05 or greater.

gene_matrix_annotated<-function(breseq,gene_list,sample_list,AF_min){
  LOR_gene_matrix<-matrix(ncol=length(gene_list),nrow=length(sample_list))
  colnames(LOR_gene_matrix)=gene_list
  rownames(LOR_gene_matrix)=sample_list
  
  for (gene in gene_list){
    for (sample in sample_list){
      if (gene%in%breseq[breseq$Sample==sample & breseq$freq>=AF_min,"gene_2"]){num_sites<-length(unique(breseq[breseq$Sample==sample & breseq$gene_2==gene & !is.na(breseq$gene_2),"position"])); LOR_gene_matrix[sample,gene]=num_sites}
      else{LOR_gene_matrix[sample,gene]=0}
    }
  }
  return(LOR_gene_matrix)
}

## Read and filter breseq output
### read files and form dataframe
breseq_annotated<-data.frame()
filenames=list.files("~/Desktop/Koskella lab/Loss of resistance/Github/Breseq")
for (file in filenames){
  ann<-data.frame(read_excel(paste("~/Desktop/Koskella lab/Loss of resistance/Github/Breseq/",file,sep="")))
  ann$Line<-unlist(strsplit(file,split="_"))[1]
  ann$Passage<-unlist(strsplit(file,split="_"))[2]
  ann$Sample<-paste(ann$Line,ann$Passage,sep="_")
  if (ncol(ann)>10){ann<-ann[,-c(8:12)]} # some files have an extra column of notes for annotations without evidence (*) resulting in blanks for cols 8-12
  breseq_annotated<-rbind(breseq_annotated,ann)
}
dim(breseq_annotated)

## Remove unassigned evidence
breseq_annotated<-breseq_annotated[!breseq_annotated$evidence%in%c(NA,"Unassigned missing coverage evidence","*","Unassigned new junction evidence","?"),]
dim(breseq_annotated)

## Remove any sites that differ from reference in ANCDC3000
ANCDC3000_sites<-breseq_annotated[breseq_annotated$Line=="ANCDC3000","position"]
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%ANCDC3000_sites,]
dim(breseq_annotated)

## Remove sites that are fixed in every single line
tab<-data.frame(table(breseq_annotated[breseq_annotated$freq==1,"Line"],breseq_annotated[breseq_annotated$freq==1,"position"])) #table of fixed sites
tab<-tab[tab$Freq>0,] # don't consider line-mutation combinations that never occurred
mut_tab<-data.frame(table(tab$Var2)) # total up number of lines in which mutation ever occurred
sites_to_remove<-as.character(mut_tab[mut_tab$Freq==length(unique(breseq_annotated$Line)),"Var1"])
breseq_annotated<-breseq_annotated[!breseq_annotated$position%in%sites_to_remove,]
dim(breseq_annotated)

## Reformat frameshift mutations to be consistent with point mutations
for (i in seq(1,nrow(breseq_annotated))){
  pos<-breseq_annotated[i,"position"]
  # find mutations with a ":"
  if (grepl(":",pos,fixed=TRUE)){
    pos<-unlist(strsplit(pos,split=":"))[1]
    pos<-paste(unlist(strsplit(pos,split=",")),collapse="")
    breseq_annotated[i,"position"]<-pos
  }
}
breseq_annotated$position<-as.numeric(breseq_annotated$position)

## Reformat gene column to remove direction
for (i in seq(1,nrow(breseq_annotated))){
  breseq_annotated[i,"gene_2"]<-substr(breseq_annotated[i,"gene"],1,nchar(breseq_annotated[i,"gene"])-2)
}

# Scan a 50 bp sliding window within each line/timepoint. If a mutation has more than N neighbors (including repeat calls at the same position), remove it.
breseq_annotated_filtered<-data.frame()
N=4

for (sample in unique(breseq_annotated$Sample)){
  pre_filter_POS<-sort(as.numeric(breseq_annotated[breseq_annotated$Sample==sample,"position"]))
  sites_to_remove<-c()
  for (pos in pre_filter_POS){
    # sliding window
    min=pos-50
    while (min<pos){
      max<-min+50
      if (length(pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max])>N)   {sites_to_remove<-unique(c(sites_to_remove,pre_filter_POS[pre_filter_POS>=min & pre_filter_POS<max]))} # if a window contains more than 3 mutations, remove all of them
      min<-min+1
    }
  }
  
  ### filter sites
  post_filter<-breseq_annotated[breseq_annotated$Sample==sample & !breseq_annotated$position%in%sites_to_remove,]
  breseq_annotated_filtered<-rbind(breseq_annotated_filtered,post_filter)
}
dim(breseq_annotated_filtered)

## Summarize annotations as NS, S, pseudogene, or intergenic
for (i in seq(1,nrow(breseq_annotated_filtered))){
  anno<-substr(breseq_annotated_filtered[i,"annotation"],1,5)
  if (anno=="pseud"){breseq_annotated_filtered[i,"Type"]="Pseudogene"}
  else if (anno=="inter"){breseq_annotated_filtered[i,"Type"]="Intergenic"}
  else {
    AA1<-substr(anno,1,1)
    AA2<-substr(anno,5,5)
    if (AA1==AA2){breseq_annotated_filtered[i,"Type"]="S"}
    else {breseq_annotated_filtered[i,"Type"]="NS"}
  }
}

## Count number of fixed sites and polymorphisms in each population at passage 12

### Fixed mutations (all types)
### Exclude resistance mutations from count (because they were acquired before experimental evolution)
res_mutations<-unique(costs_of_res$Position)
a<-aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"freq"],by=list(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"Line"]),FUN=function(x){length(x[x==1])})
mean(a$x)
sd(a$x)
t.test(x~(substr(a$Group.1,1,3)=="ANC"),a,var.equal=T)

### Polymorphisms (all types)
p<-aggregate(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"freq"],by=list(breseq_annotated_filtered[breseq_annotated_filtered$Passage=="P12" & breseq_annotated_filtered$Type%in%c("NS","S") & !breseq_annotated_filtered$position%in%c(res_mutations),"Line"]),FUN=function(x){length(x[x<1])})
mean(p$x)
sd(p$x)
t.test(x~(substr(p$Group.1,1,3)=="ANC"),p,var.equal=T)

### Fixed mutations in mutator populations (within genes)
mut<-c("FMS13","MS1","MS10","MS15","QAC4")
p[p$Group.1%in%mut,"Mutator"]<-"Y"
p[!p$Group.1%in%mut,"Mutator"]<-"N"
t.test(x~Mutator,p,alternative="less",var.equal=T)

## Figure 4: Genomic parallelism across populations
### Generate names of populations initiated with phage-resistant bacteria only (exclude ancestral controls)
n<-unique(breseq_annotated_filtered[substr(breseq_annotated_filtered$Line,1,3)!="ANC","Line"])

### Form a data frame of non-synonymous mutations not fixed at Psg 0 only
breseq_NS_noP0<-data.frame(matrix(nrow=0,ncol=ncol(breseq_annotated_filtered)))
for (line in n){
  P0_mutations<-breseq_annotated_filtered[breseq_annotated_filtered$Line==line & breseq_annotated_filtered$Passage=="P0" & breseq_annotated_filtered$freq==1,"position"]
  ###rows corresponding to P2/P12 mutations not fixed in P0 of this sample
  sample_data<-breseq_annotated_filtered[breseq_annotated_filtered$Line==line & breseq_annotated_filtered$Passage!="P0" & !breseq_annotated_filtered$position%in%P0_mutations & breseq_annotated_filtered$Type=="NS",]
  
  breseq_NS_noP0<-rbind(breseq_NS_noP0,sample_data)
}

### Make a gene matrix using only mutations not fixed in P0
gene_list<-sort(unique(breseq_NS_noP0$"gene_2"))
sample_list<-unique(breseq_NS_noP0$Sample)
LOR_gene_matrix_NS_noP0<-gene_matrix_annotated(breseq_NS_noP0,gene_list,sample_list,0)

### Calculate Sørensen-Dice similarity coefficients for every pair of samples
all_pairs<-combn(sample_list,2)
sim_coef_NS_noP0<-data.frame(matrix(nrow=0,ncol=3))

for (i in seq(1,ncol(all_pairs))){
  sample1<-all_pairs[1,i]
  sample2<-all_pairs[2,i]
  
  ### Names of genes with any mutations from ancestral genotype (any number of independent mutations per gene)
  sample1_genes<-names(LOR_gene_matrix_NS_noP0[sample1,apply(data.frame(LOR_gene_matrix_NS_noP0[sample1,]),1,function(x){x>0})])
  sample2_genes<-names(LOR_gene_matrix_NS_noP0[sample2,apply(data.frame(LOR_gene_matrix_NS_noP0[sample2,]),1,function(x){x>0})])
  
  ### Calculate the extent of overlap
  sim_coef<-(2*length(intersect(sample1_genes,sample2_genes))) / (length(sample1_genes) + length(sample2_genes))
  sim_coef_NS_noP0<-rbind(sim_coef_NS_noP0,c(sample1,sample2,sim_coef))
}
colnames(sim_coef_NS_noP0)=c("sample1","sample2","sim_coef")
sim_coef_NS_noP0$sim_coef<-as.numeric(sim_coef_NS_noP0$sim_coef)

### Annotate the resistance genes of the populations
for (i in seq(1,nrow(sim_coef_NS_noP0))){
  sim_coef_NS_noP0[i,"line1"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample1"],split="_"))[1]
  sim_coef_NS_noP0[i,"passage1"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample1"],split="_"))[2]
  sim_coef_NS_noP0[i,"line2"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample2"],split="_"))[1]
  sim_coef_NS_noP0[i,"passage2"]<-unlist(strsplit(sim_coef_NS_noP0[i,"sample2"],split="_"))[2]
  
  sim_coef_NS_noP0[i,"gene1"]=as.character(unique(costs_of_res[costs_of_res$Population==sim_coef_NS_noP0[i,"line1"],"Gene"]))
  sim_coef_NS_noP0[i,"gene2"]=as.character(unique(costs_of_res[costs_of_res$Population==sim_coef_NS_noP0[i,"line2"],"Gene"]))
  }

# Annotate whether pairs have the same or different resistance genes
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1==sim_coef_NS_noP0$gene2,"same_gene"]="same"
sim_coef_NS_noP0[!is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2) & sim_coef_NS_noP0$gene1!=sim_coef_NS_noP0$gene2,"same_gene"]="different"

## Randomization test (Day 6)
### Construct test statistic
tmp<-sim_coef_NS_noP0[sim_coef_NS_noP0$passage1=="P2" & sim_coef_NS_noP0$passage2=="P2" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2),]
test_same_P2<-tmp[tmp$same_gene=="same","sim_coef"]
test_diff_P2<-tmp[tmp$same_gene=="different","sim_coef"]
test_gene_P2<-mean(test_same_P2)-mean(test_diff_P2)

### Randomization & resampling
perm_gene_P2<-c()
set.seed(123)
for (i in seq(1,10000)){
  tmp$same_gene<-sample(tmp$same_gene,replace=F)
  perm_same_P2<-mean(tmp[tmp$same_gene=="same","sim_coef"])
  perm_diff_P2<-mean(tmp[tmp$same_gene=="different","sim_coef"])
  perm_gene_P2<-c(perm_gene_P2,perm_same_P2-perm_diff_P2)
}

### p-value (proportion of permutations more extreme than observed value)
print(length(perm_gene_P2[perm_gene_P2>=test_gene_P2])/length(perm_gene_P2))

## Randomization test (Day 36)
### Construct test statistic
tmp<-sim_coef_NS_noP0[sim_coef_NS_noP0$passage1=="P12" & sim_coef_NS_noP0$passage2=="P12" & !is.na(sim_coef_NS_noP0$gene1) & !is.na(sim_coef_NS_noP0$gene2),]
test_same_P12<-tmp[tmp$same_gene=="same","sim_coef"]
test_diff_P12<-tmp[tmp$same_gene=="different","sim_coef"]
test_gene_P12<-mean(test_same_P12)-mean(test_diff_P12)

### Randomization & resampling
perm_gene_P12<-c()
set.seed(123)
for (i in seq(1,10000)){
  tmp$same_gene<-sample(tmp$same_gene,replace=F)
  perm_same_P12<-mean(tmp[tmp$same_gene=="same","sim_coef"])
  perm_diff_P12<-mean(tmp[tmp$same_gene=="different","sim_coef"])
  perm_gene_P12<-c(perm_gene_P12,perm_same_P12-perm_diff_P12)
}
### p-value (proportion of permutations more extreme than observed value)
print(length(perm_gene_P12[perm_gene_P12>=test_gene_P12])/length(perm_gene_P12))

## Code for Figure 4A
test_gene_df<-data.frame(c(test_same_P2,test_diff_P2),c(rep("Same\nresistance\ngene",length(test_same_P2)),rep("Different\nresistance\ngene",length(test_diff_P2))))
colnames(test_gene_df)=c("sim_coef","label")
test_gene_df$label<-factor(test_gene_df$label,levels=c("Same\nresistance\ngene","Different\nresistance\ngene"))
ggplot(test_gene_df,aes(label,sim_coef))+geom_boxplot(outlier.shape=NA,width=0.4,size=0.8)+geom_jitter(width=0.05,size=2,alpha=0.2)+theme_classic(base_size=18)+xlab("")+ylab("Overlap in acquired mutations\n(Pairwise Sørensen-Dice similarity)")+theme(panel.spacing = unit(3, "lines"),strip.text=element_text(size=18,face="bold"))+ylim(0.1,0.55)

## Identify which genes were mutated in which populations
top_genes_in_study<-data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","description"],breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","Line"]))
colnames(top_genes_in_study)=c("description","Line","hits")
top_genes<-names(tail(sort(table(top_genes_in_study[top_genes_in_study$hits>0,"description"])),20))

### Construct a data frame of mutations that appeared by Day 6, excluding phage resistance mutations or any other mutations fixed before the start of the experiment
breseq_gene_analysis<-data.frame(table(breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","description"],breseq_NS_noP0[breseq_NS_noP0$Passage=="P2","Line"]))
colnames(breseq_gene_analysis)=c("gene","Line","hits")

### Annotate each population according to its resistance gene
breseq_gene_analysis[breseq_gene_analysis$Line%in%rfbA,"res_gene"]<-"rfbA"
breseq_gene_analysis[breseq_gene_analysis$Line%in%glycosyl,"res_gene"]<-"glycosyl"
breseq_gene_analysis[breseq_gene_analysis$Line%in%glycoside,"res_gene"]<-"glycoside"
breseq_gene_analysis<-breseq_gene_analysis[!is.na(breseq_gene_analysis$res_gene),]

### Calculate percent of populations with each resistance gene (rfbA, PSPTO_4988, or PSPTO_4991) that mutated each of the other gene
breseq_gene_props<-data.frame(matrix(nrow=0,ncol=3))
l<-list(rfbA,glycosyl,glycoside)
names(l)=c("rfbA","glycosyl","glycoside")
for (gene in top_genes){
  for (j in seq(1,3)){
    res_gene<-l[[j]]
    prop<-nrow(breseq_gene_analysis[breseq_gene_analysis$Line%in%res_gene & breseq_gene_analysis$gene==gene & breseq_gene_analysis$hits>0,])/nrow(breseq_gene_analysis[breseq_gene_analysis$Line%in%res_gene & breseq_gene_analysis$gene==gene,])
    breseq_gene_props<-rbind(breseq_gene_props,c(gene,names(l)[[j]],prop))
  }
}
colnames(breseq_gene_props)=c("top_gene","res_gene","percent_pops")
breseq_gene_props$percent_pops<-as.numeric(breseq_gene_props$percent_pops)
breseq_gene_props[breseq_gene_props$percent_pops==0,"percent_pops"]<-10^-4

# Fig 4C (9x6)
breseq_gene_props$top_gene<-factor(breseq_gene_props$top_gene,levels=rev(unique(breseq_gene_props$top_gene)))
ggplot(breseq_gene_props,aes(top_gene,res_gene,fill=percent_pops*100))+geom_tile(color="grey80")+theme_classic(base_size=14)+theme(axis.text.x=element_text(angle=90,hjust=1))+scale_fill_viridis()+ylab("")+xlab("")

