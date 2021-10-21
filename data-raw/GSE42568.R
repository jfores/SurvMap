## code to prepare `GSE42568` dataset goes here

rm(list = ls())
options(stringsAsFactors = F)

#######################################
#Bash code to get the dataset from GEO#
#######################################

#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/suppl/GSE42568_RAW.tar
#wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE42nnn/GSE42568/soft/GSE42568_family.soft.gz
#tar -xf GSE42568_RAW.tar
#gunzip *.CEL.gz

###################
#Loading libraries#
###################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligo",force = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("oligoClasses",force = TRUE)

install.packages("ff")
install.packages("/home/jj/Downloads/ff_2.2-14.tar.gz", repos = NULL, type="source")
install.packages("/home/jj/Downloads/bit_1.1-15.2.tar.gz", repos = NULL, type="source")


library(affy)
library(oligoClasses)
library(oligo)
library(GEOquery)
library(DupChecker)
library(arrayQualityMetrics)
library(vctrs)
library(frma)
library(a4Preproc)
library(genefu)
library(plyr)
library(preprocessCore)

#####################################
#Loading and reformating sample data#
#####################################

GEO_dat  <- getGEO(filename= "/home/jj/Desktop/temp/GSE42568/Data/GSE42568_family.soft.gz")
platforms <- unlist(lapply(GPLList(GEO_dat),function(x) Meta(x)$geo_accession))
pData_GEO_dat <- get_pData_table(GEO_dat,platforms[[1]])
rownames(pData_GEO_dat) <- pData_GEO_dat$geo_accession
days_to_months <- 30.4375
pCh_Status <- pData_GEO_dat$source_name_ch1
pCh_Status[grepl("cancer",pCh_Status)] <- "T"
pCh_Status[grepl("normal",pCh_Status)] <- "NT"
pCh_Age <- pData_GEO_dat$pCh_age
pCh_Node <- pData_GEO_dat$pCh_lymph_node_status
pCh_Size <-  pData_GEO_dat$pCh_size
pCh_ER_Status <- pData_GEO_dat$pCh_er_status
pCh_Grade <- pData_GEO_dat$pCh_grade
pCh_DFS_T <- pData_GEO_dat$pCh_relapse_free_survival_time_days
pCh_DFS_T <- as.numeric(pCh_DFS_T)/days_to_months
pCh_DFS_E <- pData_GEO_dat$pCh_relapse_free_survival_event
pCh_Death_T <- pData_GEO_dat$pCh_overall_survival_time_days
pCh_Death_T <- as.numeric(pCh_Death_T)/days_to_months
pCh_Death_E <- pData_GEO_dat$pCh_overall_survival_event

vars <- ls()[grepl("pCh_",ls())]
vars <- vars[vars != "pCh_to_Num_and_Fac"]
vars <- vars[vars != "apply_NA_To_Undet_To_pCh_Cols"]

eval(parse(text = paste("pData_GEO_ft <- data.frame(",paste(vars,collapse = ","),")",sep="")))
rownames(pData_GEO_ft) <- rownames(pData_GEO_dat)


#Removing cel files generated with the B platform.

files_temp <- dir("/home/jj/Desktop/temp/GSE42568/Data/",full.names = T)
file.remove(files_temp[grepl("soft.gz",files_temp)])
file.remove(files_temp[grepl("_RAW.tar",files_temp)])
file.remove(files_temp[grepl("CEL.md5",files_temp)])

#Obtaining md5 for each .CEL file.

datafile <- buildFileTable(rootDir="/home/jj/Desktop/temp/GSE42568/Data/")
result <-validateFile(datafile)
result <- result$filetable
rownames(result) <- gsub("_.*","",gsub(".CEL","",gsub(".*/Data/","",result$file)))
result$pCh_samples_to_df <- rownames(result)
result <- result[,-1]
result <- result[,c(3,1,2)]
colnames(result) <- c("pCh_samples_to_df","pCh_directory_to_df","pCh_md5_to_df")
df_out <- result

#Merging

pData_GEO_ft <- merge(pData_GEO_ft,df_out,by = "row.names")
colnames(pData_GEO_ft)[1] <- "pCh_Sample_Name"
rownames(pData_GEO_ft) <- pData_GEO_ft[,1]

directory_files_temp <- dir("/home/jj/Desktop/temp/GSE42568/Data/",full.names = T)
directory_files_temp_to_delete <- directory_files_temp[grepl(".CEL.md5$",directory_files_temp)]
for(i in 1:length(directory_files_temp_to_delete)){
  file.remove(directory_files_temp_to_delete[i])
}

#Reading files.

directory_files_temp <- dir("/home/jj/Desktop/temp/GSE42568/Data/",full.names = T)
eData_GEO <- ReadAffy(filenames=directory_files_temp)
cols_eDat <- colnames(eData_GEO)
cols_eDat <- gsub("_.*","",gsub(".CEL","",cols_eDat))
pData_GEO_ft <- pData_GEO_ft[cols_eDat,]

#Generating expressionset

GEO_Eset_Raw <- affy::rma(eData_GEO,background = TRUE,normalize = TRUE)
help(rma)
colnames(GEO_Eset_Raw)
colnames(GEO_Eset_Raw) <- gsub("_.*","",gsub(".CEL","",colnames(GEO_Eset_Raw)))
pData(GEO_Eset_Raw) <- pData_GEO_ft[colnames(GEO_Eset_Raw),]
exprs(GEO_Eset_Raw)
BiocManager::install("preprocessCore", configure.args="--disable-threading")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages("ropenblas")
BiocManager::install("oligo",force = T)

###########
#Functions#
###########


unfold_charact <- function(GSM){
  GSM_meta <- Meta(GSM)
  characteristics <- GSM_meta$characteristics_ch1
  for(i in 1:length(characteristics)){
    if(length(strsplit(characteristics[i],":")[[1]]) > 1){
      GSM_meta[[paste("pCh_",gsub(" ","_",trimws(strsplit(characteristics[i],":")[[1]][1])),sep="")]] <- gsub(" ","_",trimws(strsplit(characteristics[i],":")[[1]][2]))
    }
  }
  GSM_meta[["characteristics_ch1"]] <- NULL
  for(i in 1:length(GSM_meta)){
    if(length(GSM_meta[[i]])>1){
      GSM_meta[[i]] <- paste(GSM_meta[[i]],collapse="/")
    }
  }
  return(GSM_meta)
}

#getting pheno data from soft

get_pData_table <- function(GSE,GPL){
  counter_temp <- 0
  for(i in 1:length(GSMList(GSE))){
    if(Meta(GSMList(GSE)[[i]])$platform_id == GPL){
      if(counter_temp == 0){
        df_sal <- data.frame(unfold_charact(GSMList(GSE)[[i]]))
      }else{
        df_sal <- rbind.fill(df_sal,data.frame(unfold_charact(GSMList(GSE)[[i]])))
      }
      counter_temp = counter_temp + 1
    }
  }
  return(df_sal)
}


usethis::use_data(GSE42568, overwrite = TRUE)

