library("RcppArmadillo")
library("Rcpp")
library("RcppParallel")
library("devtools")
devtools::install_github("stop-pre16/mcmseq", build_vignettes = T)
library("mcmseq")

gene2018 <- "light_dark_rna_fpkms_CPM_DEG_18june19_count.csv" #file containing the read count matrix


cts <- as.matrix(read.csv(gene2018,row.names="gene", header = TRUE))

cts=cts[,-c(1:4)]#retain only counts, remove all other gene info
dim(cts)
cts_num=t(apply(cts, 1,function(x) as.numeric(as.character(x))))
colnames(cts_num)=colnames(cts)

dim(cts_num)

ind=c()
for (i in 1:nrow(cts_num)){
  if(all(cts_num[i,] == 0)){ind=c(ind,i)}
}
countsNonZero <- cts_num[-ind,]

#How many rows are left?
nrow(countsNonZero)
colnames(countsNonZero)
countsNonZero <- countsNonZero[,order(colnames(countsNonZero))]


# Generating the coldata table -------------------------------------------------------------------

coldata=data.frame(matrix(colnames(countsNonZero)))
coldata['rep']=as.factor(gsub('.*_.*_R','',as.character(coldata[,1])))
coldata['treatment']= as.factor(gsub('.*_(.*)_.*','\\1',as.character(coldata[,1])))
coldata['genotype']=as.factor(gsub('_.*_.*','',as.character(coldata[,1])))
coldata["geno_trt"] <- paste(coldata$genotype, coldata$treatment, sep = "_")
coldata["geno_rep"] <- paste(coldata$genotype, coldata$rep, sep = "_")
coldata <- coldata[,-1]
rownames(coldata)


f <- ~ rep + genotype:treatment

x <- model.matrix(f, data = coldata)
colnames(x)

contrast <-rbind(c(0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0,0),
                 c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0),
                 c(0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0,0),
                 c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0),
                 c(0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0),
                 c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0),
                 c(0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0),
                 c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0),
                 c(0,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0,0),
                 c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0),
                 c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1,0),
                 c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0),
                 c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,-1,0,0,0,0,0,0,0),
                 c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,1),
                 c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,-1))




rownames(contrast) <- c("B73:c-d", "B73:l-d", "B73:c-l","B97:c-d", "B97:l-d", "B97:c-l","Ki11:c-d", "Ki11:l-d", "Ki11:c-l","M37W:c-d", "M37W:l-d", "M37W:c-l","MS71:c-d", "MS71:l-d", "MS71:c-l","NC358:c-d", "NC358:l-d", "NC358:c-l","OH7B:c-d", "OH7B:l-d", "OH7B:c-l")

ncol(contrast)

fit.default <- mcmseq.fit(counts=countsNonZero, 
                          fixed_effects = f,
                          sample_data = coldata,
                          random_intercept = 'geno_rep', 
                          gene_names = rownames(countsNonZero),
                          contrast_mat = contrast,
                          contrast_names = NULL,
                          n_it = 50000,
                          prop_burn_in = 0.1
)

B73_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                         summarizeWhat="contrast",  
                         which = "B73:c-d",
                         order_by = "BH_adjusted_pval", 
                         decreasing = FALSE,
                         filter_by = c("BH_adjusted_pval"),
                         filter_val = c(0.05),
                         
                         log2=TRUE
)

B73_CD$contrast <- rep("B73_CD",lengths(B73_CD)[1])

B73_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                         summarizeWhat="contrast",  
                         which = "B73:c-d",
                         order_by = "BH_adjusted_pval", 
                         decreasing = FALSE,
                         filter_by = c("BH_adjusted_pval"),
                         filter_val = c(0.05),
                         
                         log2=TRUE
)

B73_LD$contrast <- rep("B73_LD",lengths(B73_LD)[1])

B73_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                         summarizeWhat="contrast",  
                         which = "B73:c-d",
                         order_by = "BH_adjusted_pval", 
                         decreasing = FALSE,
                         filter_by = c("BH_adjusted_pval"),
                         filter_val = c(0.05),
                         
                         log2=TRUE
)
B73_CL$contrast <- rep("B73_CL",lengths(B73_CL)[1])


B97_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                         summarizeWhat="contrast",  
                         which = "B97:c-d",
                         order_by = "BH_adjusted_pval", 
                         decreasing = FALSE,
                         filter_by = c("BH_adjusted_pval"),
                         filter_val = c(0.05),
                         
                         log2=TRUE
)

B97_CD$contrast <- rep("B97_CD",lengths(B97_CD)[1])


B97_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                         summarizeWhat="contrast",  
                         which = "B97:c-d",
                         order_by = "BH_adjusted_pval", 
                         decreasing = FALSE,
                         filter_by = c("BH_adjusted_pval"),
                         filter_val = c(0.05),
                         
                         log2=TRUE
)

B97_LD$contrast <- rep("B97_LD",lengths(B97_LD)[1])


B97_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                         summarizeWhat="contrast",  
                         which = "B97:c-d",
                         order_by = "BH_adjusted_pval", 
                         decreasing = FALSE,
                         filter_by = c("BH_adjusted_pval"),
                         filter_val = c(0.05),
                         
                         log2=TRUE
)

B97_CL$contrast <- rep("B97_CL",lengths(B97_CL)[1])


Ki11_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "Ki11:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)
Ki11_CD$contrast <- rep("Ki11_CD",lengths(Ki11_CD)[1])

Ki11_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "Ki11:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)
Ki11_LD$contrast <- rep("Ki11_LD",lengths(Ki11_LD)[1])


Ki11_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "Ki11:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)


Ki11_CL$contrast <- rep("Ki11_CL",lengths(Ki11_CL)[1])


M37W_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "M37W:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)

M37W_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "M37W:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c( 0.05),
                          
                          log2=TRUE
)
M37W_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "M37W:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)
MS71_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "MS71:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)

MS71_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "MS71:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)
MS71_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "MS71:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)

NC358_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                           summarizeWhat="contrast",  
                           which = "NC358:c-d",
                           order_by = "BH_adjusted_pval", 
                           decreasing = FALSE,
                           filter_by = c("BH_adjusted_pval"),
                           filter_val = c(0.05),
                           
                           log2=TRUE
)

NC358_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                           summarizeWhat="contrast",  
                           which = "NC358:c-d",
                           order_by = "BH_adjusted_pval", 
                           decreasing = FALSE,
                           filter_by = c("BH_adjusted_pval"),
                           filter_val = c(0.05),
                           
                           log2=TRUE
)
NC358_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                           summarizeWhat="contrast",  
                           which = "NC358:c-d",
                           order_by = "BH_adjusted_pval", 
                           decreasing = FALSE,
                           filter_by = c("BH_adjusted_pval"),
                           filter_val = c(0.05),
                           
                           log2=TRUE
)

OH7B_CD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "OH7B:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)

OH7B_LD <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "OH7B:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)
OH7B_CL <- mcmseq.summary(mcmseqModel = fit.default, 
                          summarizeWhat="contrast",  
                          which = "OH7B:c-d",
                          order_by = "BH_adjusted_pval", 
                          decreasing = FALSE,
                          filter_by = c("BH_adjusted_pval"),
                          filter_val = c(0.05),
                          
                          log2=TRUE
)


OH7B_CD$contrast <- rep("OH7B_CD",lengths(OH7B_CD)[1])
OH7B_LD$contrast <- rep("OH7B_LD",lengths(OH7B_LD)[1])
OH7B_CL$contrast <- rep("OH7B_CL",lengths(OH7B_CL)[1])

M37W_CD$contrast <- rep("M37W_CD",lengths(M37W_CD)[1])
M37W_LD$contrast <- rep("M37W_LD",lengths(M37W_LD)[1])
M37W_CL$contrast <- rep("M37W_CL",lengths(M37W_CL)[1])

MS71_CD$contrast <- rep("MS71_CD",lengths(MS71_CD)[1])
MS71_LD$contrast <- rep("MS71_LD",lengths(MS71_LD)[1])
MS71_CL$contrast <- rep("MS71_CL",lengths(MS71_CL)[1])

NC358_CD$contrast <- rep("NC358_CD",lengths(NC358_CD)[1])
NC358_LD$contrast <- rep("NC358_LD",lengths(NC358_LD)[1])
NC358_CL$contrast <- rep("NC358_CL",lengths(NC358_CL)[1])



df <- rbind(B73_CD, B73_CL, B73_LD,B97_CD, B97_CL,B97_LD, Ki11_CD, Ki11_CL, Ki11_LD, M37W_CD, M37W_CL, M37W_LD, MS71_CD, MS71_CL, MS71_LD, NC358_CD, NC358_CL, NC358_LD, OH7B_CD,OH7B_CL, OH7B_LD)

write.csv(df, "MCMseq_output.csv", row.names =  F)

