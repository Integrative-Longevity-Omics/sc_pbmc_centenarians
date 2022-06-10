#!/usr/bin/env Rscript

#arguments from commandline
args = commandArgs(trailingOnly=TRUE)
n.ct <- args[1]
n.cores <- args[2]
n.genes <- args[3]

n.ct <- as.numeric(n.ct)
n.cores <- as.numeric(n.cores)
n.genes <- as.numeric(n.genes)

print(paste0("cell.type=",n.ct))
print(paste("number of cores =", n.cores))
print(paste("subset of genes =", n.genes))

#load libraries and additional functions
library(rjags)
library(coda)
library(dplyr)
library(Seurat)
library(parallel)
library(scater)
library(reshape2)
library(stringr)

#Load data
work.dir <- "/restricted/projectnb/uh2-sebas/analysis/scCentenarians/"

#source model
source(paste0(work.dir, "Differential_Analysis/Bayes_Model/LifeSpan_Model/jags.model.R"))

pbmc.combined <- readRDS(file = paste0(work.dir,"Integrated_Data/pbmc.combined.rds"))


#keep only immune cells
Idents(pbmc.combined) <- "ct.consensus"
pbmc.combined <- subset(pbmc.combined, idents = c("Plasma", "CD4 T Cytotoxic", "T gamma delta","pDC"), invert = TRUE)
pbmc.combined$ct.consensus <- factor(pbmc.combined$ct.consensus, levels = c("CD4 T Naive", "CD4 T Memory", "NK", "M14", "BC Memory","BC Naive","CD8 T","M16", "mDC"))
cell.types <- levels(pbmc.combined$ct.consensus)
print(cell.types)

#filtering genes to include genes expressed in over 50% of cells of the least frequent cell type
min.count <- sum(pbmc.combined@meta.data$ct.consensus == "mDC")/2

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- pbmc.combined[['RNA']]@counts > 0

# Sums all TRUE values and returns TRUE if more than min.count TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= min.count
length(keep_genes)

#filter genes based on threshold
pbmc.combined[["RNA"]]@data <- pbmc.combined[["RNA"]]@data[keep_genes,]
dim(pbmc.combined)

#subset data on cell type
Idents(pbmc.combined) <- "ct.consensus"
pbmc.set <- subset(pbmc.combined, idents = c(cell.types[n.ct])) 
pbmc.set$ct.consensus <- droplevels(pbmc.set$ct.consensus) 
rm(pbmc.combined)
print(paste("cell type of interest:", levels(pbmc.set$ct.consensus)))

#dummy variable for age groups
#four age groups
pbmc.set@meta.data$age.quartiles <- factor(pbmc.set@meta.data$age.quartiles, levels = c("20-39", "40-59", "60-89","Extreme Longevity"))
library(fastDummies)
age.dummy <- dummy_cols(pbmc.set@meta.data$age.quartiles)
pbmc.set@meta.data$age.group.middle <- age.dummy$`.data_40-59`
pbmc.set@meta.data$age.group.old <- age.dummy$`.data_60-89`
pbmc.set@meta.data$age.group.EL <- age.dummy$`.data_Extreme Longevity`

#dummy variable for ethnicity
pbmc.set@meta.data$ethnic[pbmc.set@meta.data$batch == "NATGEN"] <- 0
pbmc.set@meta.data$ethnic[pbmc.set@meta.data$batch == "COHORT.B1"] <- 0
pbmc.set@meta.data$ethnic[pbmc.set@meta.data$batch == "COHORT.B2"] <- 0
pbmc.set@meta.data$ethnic[pbmc.set@meta.data$batch == "PNAS"] <- 1

#dummy variable for batch
pbmc.set@meta.data$batch[pbmc.set@meta.data$batch == "NATGEN"] <- 1
pbmc.set@meta.data$batch[pbmc.set@meta.data$batch == "COHORT.B1"] <- 2
pbmc.set@meta.data$batch[pbmc.set@meta.data$batch == "COHORT.B2"] <- 3
pbmc.set@meta.data$batch[pbmc.set@meta.data$batch == "PNAS"] <- 4
pbmc.set@meta.data$batch <- as.factor(pbmc.set@meta.data$batch)

#dummy variable for sex
pbmc.set@meta.data$sex <- as.factor(pbmc.set@meta.data$sex)
pbmc.set@meta.data$sex <- relevel(pbmc.set@meta.data$sex, "Male")
pbmc.set@meta.data$sex <- as.numeric(pbmc.set@meta.data$sex == "Female")    

#dummy variable for sample ID
pbmc.set@meta.data$sample.ID <- as.integer(factor(pbmc.set@meta.data$sample.ID))


#Ceiling
hvgenes <- rownames(pbmc.set)
SC.chunks <- split(1:length(hvgenes), ceiling(seq_along(1:length(hvgenes))/(length(hvgenes)/100)))
print(SC.chunks[[n.genes]])
chunk <- SC.chunks[[n.genes]]

hvgenes <- hvgenes[chunk]

# print("created and selected partition")
    
#run bayesian model through loop
cl <- makeCluster(n.cores)
clusterEvalQ(cl, {
	  library(rjags) 
	  library(coda)
	})
clusterExport(cl, list("pbmc.set","hvgenes", "model"))
start.time <- Sys.time()
	
jags.results <- parLapply(cl, X = 1:length(hvgenes), function(i){
	  
	  data <- list(y.gene = as.matrix(pbmc.set[['RNA']]@data[hvgenes[i],])[,1], #expression across cells
	               N.cell = ncol(pbmc.set), #total number of cells
	               N.subj = length(levels(as.factor(pbmc.set@meta.data$sample.ID))), #number of samples
	               subj = pbmc.set@meta.data$sample.ID, #sample ID for each cell i.e. factored 1,2,..n
	               age.group.middle = pbmc.set@meta.data$age.group.middle, #middle age group v. young age group i.e. binary
	               age.group.old = pbmc.set@meta.data$age.group.old, #old age group v. young age group i.e. binary
	               age.group.EL = pbmc.set@meta.data$age.group.EL, #EL group v. young age group i.e. binary
	               sex = pbmc.set@meta.data$sex, #sex for each cell i.e. binary
                 ethnic = pbmc.set@meta.data$ethnic, #ethnicity for each cell i.e. binary
	               batch = pbmc.set@meta.data$batch #batch ID for each cell i.e. factored 1,2,..n
	  )
	  set.seed(123)
	  jags <- jags.model(textConnection(model),data=data, n.adapt=1500)
	  test <- coda.samples(jags, c('b1','b2','b3','bs','tau','tau.s'), n.adapt = 1500, n.iter=10000)
	  
	  coda.summary <- summary(test[,c('b1','b2','b3')])
      coda.summary <- data.frame(as.matrix(coda.summary$statistics), as.matrix(coda.summary$quantiles[,c("2.5%","97.5%")]))
      rownames(coda.summary) <- c("Middle.v.Young", "Old.v.Young","EL.v.Young")
      
      y <- abs(geweke.diag(test)[[1]][[1]])
      if (any(y[c("b1", "b2","b3")] > 2)){
      #return 0 if have an absolute value greater than 2
      coda.summary$converged <- 0
      } else {
      #return 1 if have an absolute value between -2 and 2
      coda.summary$converged <- 1
      }
      
      #subject within variance
      subj.var <- summary(test[,c('tau','tau.s')])
      subj.var <- data.frame(as.matrix(subj.var$statistics), as.matrix(subj.var$quantiles[,c("2.5%","97.5%")]))
      coda.summary$subj.var <- (subj.var["tau","Mean"] - subj.var["tau.s", "Mean"]) / subj.var["tau.s","Mean"]
      
	  return(list(test, coda.summary))
	}
	)
	
	
end.time <- Sys.time()
cat('Time to fit model',(end.time - start.time))
stopCluster(cl)
    
names(jags.results) <- hvgenes
    
print("end model")
    
doc.name <- str_replace_all(cell.types[n.ct], fixed(" "), "_")
dir.create(paste0(work.dir,"Differential_Analysis/JAGS_Output/JAGS.Lifespan/", doc.name)) 
saveRDS(jags.results, paste0(work.dir,"Differential_Analysis/JAGS_Output/JAGS.Lifespan/",doc.name,"/jags.results.",n.genes,".rds"))
    
