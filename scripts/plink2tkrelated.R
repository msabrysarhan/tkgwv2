### Estimating relatedness from individual text-PLINKs (Parallel Version)
### Author - Daniel Fernandes
# Programs with system-wide installation required # 
# 1. plink
# R packages required # 
# 1. data.table
# 2. parallel
# 3. doParallel
# 4. foreach

# Description of arguments #
# dyads = Default TRUE will analyse all possible dyads between all individuals in current folder.
#         Otherwise, user can input 'dyads' as a tab-spaced text file with each line containing a pair to be analysed.
#         Sample1 Sample2
#         Sample1 Sample3
#         (...)
# freqFile = Allele frequencies file in PLINK format (FRQ) for the same SNP set used in bam2plink() (or for example for the 1240K dataset)
# ignoreThresh = Default 1. Threshold for the minimum number of SNPs allowed to estimate relatedness.
# verbose = Default FALSE. Use verbose mode and print extended information throughout file processing.
# threads = Number of CPU threads to use (default = 1)
#
# Input:
# - individual text-PLINK
# - pairwise FRQ
# Output:
# - pairwise relatedness coefficient result text file  

# Load required packages
library(data.table)
library(parallel)
library(doParallel)
library(foreach)

args = commandArgs(trailingOnly=TRUE)
pywd = args[length(args)]
args = args[-length(args)]

if(length(grep("bam2plink",args)) != 0) { 
  args = args[-grep("bam2plink",args)]
} else { 
  cat("\n ################################################################################")
  cat("\n ### TKGWV2 - An ancient DNA relatedness pipeline for ultra-low coverage data ###")
  cat("\n ## Version 1.0b - Released 07/2022 (Parallel Version)")
  cat("\n #")
}

cat(paste0("\n # [",Sys.time(),"] ","Running 'plink2tkrelated' on folder ", getwd()))
cat("\n\t # Text-PLINK >> Pairwise transposed text-PLINK >> Relatedness estimates ")

# Parse arguments
if(length(grep("freqFile",args)) != 0) {
  freqFile = strsplit(args[grep("freqFile",args)],"=")[[1]][2]
} else {
  cat("\n\t # ERROR: The required argument 'freqFile' could not be found\n")
  quit()
}

dyads = TRUE
if(length(grep("dyads",args)) != 0) { dyads = strsplit(args[grep("dyads",args)],"=")[[1]][2] }

ignoreThresh = 1
if(length(grep("ignoreThresh",args)) != 0) { ignoreThresh = strsplit(args[grep("ignoreThresh",args)],"=")[[1]][2] }

verbose = FALSE
if(length(grep("verbose",args)) != 0) { verbose = TRUE }

# Add thread parameter
threads = 1
if(length(grep("threads",args)) != 0) { 
  threads = as.numeric(strsplit(args[grep("threads",args)],"=")[[1]][2]) 
  # Limit threads to available cores
  available_cores <- detectCores()
  if(threads > available_cores) {
    cat(paste0("\n\t # Warning: Requested ", threads, " threads but only ", available_cores, " available. Using ", available_cores, " threads.\n"))
    threads = available_cores
  }
}

# Function to process a single dyad
process_dyad <- function(line, combos, freqFile, ignoreThresh, verbose, sampleList) {
  samp1 = combos[line,1]
  samp2 = combos[line,2]
  
  cat(paste0("\t # Estimating coefficient of relatedness Rxy for   ",samp1,"   ",samp2,"\t"))
  
  # Create temp files with unique names to avoid collisions
  temp_prefix <- paste0("temp_", Sys.getpid(), "_", line, "_")
  comm_snps_file <- paste0(temp_prefix, "commSNPs")
  samp1_short <- paste0(temp_prefix, samp1, "short")
  samp2_short <- paste0(temp_prefix, samp2, "short")
  merged_prefix <- paste0(temp_prefix, samp1, "____", samp2)
  comm_frq_file <- paste0(temp_prefix, "comm.frq")
  
  #Find common SNPs
  comm4 = paste0("grep -F -x -f ",samp1,".map ",samp2,".map | awk '{print $2}' > ", comm_snps_file)
  system(comm4)
  
  result <- data.frame(Sample1=samp1, Sample2=samp2, Used_SNPs=NA, HRC=NA, 
                      counts0=NA, counts4=NA, Relationship=NA, stringsAsFactors=FALSE)
  
  if(file.size(comm_snps_file) > 0) {
    #Use PLINK to reduce samples to common SNPs and merge
    comm5 = paste0("plink --file ",samp1," --extract ", comm_snps_file," --make-bed --out ", samp1_short," --allow-no-sex")
    system(comm5, ignore.stdout = TRUE, ignore.stderr = TRUE)
    comm6 = paste0("plink --file ",samp2," --extract ", comm_snps_file," --make-bed --out ", samp2_short," --allow-no-sex")
    system(comm6, ignore.stdout = TRUE, ignore.stderr = TRUE)
    comm7 = paste0("plink --bfile ", samp1_short," --bmerge ", samp2_short," --geno 0.1 --recode transpose --out ", merged_prefix," --allow-no-sex")
    system(comm7, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if(file.exists(paste0(merged_prefix,".tped"))) {
      comm8 = paste0("awk '{print $2}' ", merged_prefix,".tped > ", comm_snps_file, "2")
      system(comm8)
      
      #Get frequencies from .frq file
      comm9 = paste0("grep -w -F -f ", comm_snps_file, "2 ", freqFile," > ", comm_frq_file)
      system(comm9)
      
      # Check for mismatches
      comm10a = paste0("wc -l ", comm_frq_file)
      comm10b = paste0("wc -l ", merged_prefix,".tped")
      if(as.integer(strsplit(system(comm10a, intern=T)," ")[[1]][1]) != as.integer(strsplit(system(comm10b, intern=T)," ")[[1]][1])) {
        comm10c = paste0("awk '{print $2\" \"}' ", comm_frq_file," > ", comm_snps_file, "3")
        system(comm10c)
        comm10d = paste0("plink --tfile ", merged_prefix," --extract ", comm_snps_file, "3 --recode transpose --out ", merged_prefix," --allow-no-sex")
        system(comm10d, ignore.stdout = TRUE, ignore.stderr = TRUE)
      }
      
      ## Read pairwise .TPED
      pairLoadNew = tryCatch({
        read.table(paste0(merged_prefix,".tped"),header=FALSE,stringsAsFactors = FALSE, sep=" ")
      }, error = function(e) {
        NULL
      })
      
      if(!is.null(pairLoadNew) && ncol(pairLoadNew) == 8) {
        pairLoadNew = pairLoadNew[,-3]
        
        ## Read allele frequencies
        alFreqNew = tryCatch({
          read.csv(comm_frq_file, stringsAsFactors=FALSE, sep="", 
                  colClasses = c("character","character","character","character","numeric","numeric"),
                  header = F, col.names = c("CHR","SNP","A1","A2","MAF","NCHROBS"))
        }, error = function(e) {
          NULL
        })
        
        if(!is.null(alFreqNew)) {
          alFreqNew$NCHROBS=NULL
          alFreqNew[5]=lapply(alFreqNew[5],round,5)
          alFreqNew$AFa2 = as.numeric(0)
          alFreqNew$AFa2=(1-alFreqNew$MAF)
          alFreqNew[6]=lapply(alFreqNew[6],round,5)
          
          ## Add allele frequency data to main dataframe
          if((length(pairLoadNew$V2) == length(alFreqNew$SNP))) {
            pairLoadNew$A1Mi = alFreqNew$A1
            pairLoadNew$A2Ma = alFreqNew$A2
            pairLoadNew$al1freq = alFreqNew$MAF
            pairLoadNew$al2freq = alFreqNew$AFa2
          } 
          
          ## Remove SNPs with fixed alleles
          pairLoadNew = pairLoadNew[pairLoadNew$al1freq != 0,]
          ## Remove SNPs with no data on .FRQ file
          if(length(which(is.na(pairLoadNew$al1freq)) != 0)) {  
            pairLoadNew = pairLoadNew[-(which(is.na(pairLoadNew$al1freq))),]  
          }
          
          ## Make sure all variants are SNPs
          rem1 = which(nchar(alFreqNew$A1) >1)
          rem2 = which(nchar(alFreqNew$A2) >1)
          rem3 = c(rem1,rem2)
          remSnps = alFreqNew[rem3,2]
          if(length(rem1) >0) { alFreqNew = alFreqNew[-rem1,] }
          if(length(rem2) >0) { alFreqNew = alFreqNew[-rem2,] }
          rem10 = which(pairLoadNew$V2 %in% remSnps)
          if(length(rem10) >0){  pairLoadNew = pairLoadNew[-rem10,] }
          rem4 = which(nchar(pairLoadNew$V5) >1)
          rem5 = which(nchar(pairLoadNew$V6) >1)
          rem6 = which(nchar(pairLoadNew$V7) >1)
          rem7 = which(nchar(pairLoadNew$V8) >1)
          rem8 = c(rem4,rem5,rem6,rem7)
          if(length(rem8) >0){  pairLoadNew = pairLoadNew[-rem8,] }
          
          ## Test if both samples have all homozygous SNPs, and force so if FALSE
          if(all(pairLoadNew$V5 == pairLoadNew$V6) == FALSE) {
            setDT(pairLoadNew)
            pairLoadNew[, V5 := V5][runif(.N, 0, 1) < 0.5 , V5 := V6][]
            pairLoadNew[, V6 := V5][runif(.N, 0, 1) < 0.5 , V6 := V6][]
            setDF(pairLoadNew)
          }
          if(all(pairLoadNew$V7 == pairLoadNew$V8) == FALSE) {
            setDT(pairLoadNew)
            pairLoadNew[, V7 := V7][runif(.N, 0, 1) < 0.5 , V7 := V8][]
            pairLoadNew[, V8 := V7][runif(.N, 0, 1) < 0.5 , V8 := V8][]
            setDF(pairLoadNew)
          }
          
          ## Work with 5 decimal places for frequencies
          maf=format(as.numeric(pairLoadNew$al1freq),signif=5)
          pairLoadNew$al1freq=maf
          af2=format((1-as.numeric(maf)),signif=5)
          pairLoadNew$al2freq=af2
          
          if(ncol(pairLoadNew) == 11) {
            colnames(pairLoadNew) = c("chr","snp","pos","S1x","S1y","S2x","S2y","A1Mi","A2Ma","al1freq","al2freq")
            
            ## Remove triallelic SNPs
            toRemove=c()
            ## Check in individual 1
            it=1
            for(i in pairLoadNew$S1x) {
              if((i!=pairLoadNew$A1Mi[it])==TRUE && (i!=pairLoadNew$A2Ma[it])==TRUE) {
                toRemove=c(toRemove,pairLoadNew$snp[it])
              }
              it=it+1}
            ## Check in individual 2
            it=1
            for(i in pairLoadNew$S2x) {
              if((i!=pairLoadNew$A1Mi[it])==TRUE && (i!=pairLoadNew$A2Ma[it])==TRUE) {
                toRemove=c(toRemove,pairLoadNew$snp[it])
              }
              it=it+1}
            ## Remove those SNPs
            toRemove = unique(toRemove)
            poses = which(pairLoadNew$snp %in% toRemove)
            if(!is.null(toRemove)) {
              pairLoadNew=pairLoadNew[-poses,]
            }
            
            result$Used_SNPs = length(pairLoadNew$snp)
            
            if(length(pairLoadNew$chr) >= ignoreThresh) {
              pairLoadNewRxy = pairLoadNew
              ## Remove duplicated allele
              pairLoadNewRxy = pairLoadNewRxy[,-5]
              pairLoadNewRxy = pairLoadNewRxy[,-6]
              ## Calculate components
              pairLoadNewRxy["PaXy"] = ifelse(as.character(pairLoadNewRxy$S1x) == pairLoadNewRxy$A1Mi, pairLoadNewRxy$al1freq, pairLoadNewRxy$al2freq)
              pairLoadNewRxy["PcXy"] = ifelse(as.character(pairLoadNewRxy$S2x) == pairLoadNewRxy$A1Mi, pairLoadNewRxy$al1freq, pairLoadNewRxy$al2freq)
              pairLoadNewRxy["IacADbcBD"] = ifelse(as.character(pairLoadNewRxy$S1x) == as.character(pairLoadNewRxy$S2x), 1, 0)*4
              pairLoadNewRxy$rxyL = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PaXy)*2)) / (1+1-(as.numeric(pairLoadNewRxy$PaXy)*2))
              pairLoadNewRxy$ryxL = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PcXy)*2)) / (1+1-(as.numeric(pairLoadNewRxy$PcXy)*2))
              rxyL = sum(pairLoadNewRxy$rxyL)
              ryxL = sum(pairLoadNewRxy$ryxL)
              pairLoadNewRxy$Av = (pairLoadNewRxy$rxyL + pairLoadNewRxy$ryxL) /2
              pairLoadNewRxy$SrQ1 = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PaXy)*2))
              pairLoadNewRxy$SrQ2 = ((0.5*(as.numeric(pairLoadNewRxy$IacADbcBD)))-(as.numeric(pairLoadNewRxy$PcXy)*2))
              pairLoadNewRxy$SrQ1MinusSum = pairLoadNewRxy$SrQ1 - sum(pairLoadNewRxy$SrQ1)
              pairLoadNewRxy$SrQ2MinusSum = pairLoadNewRxy$SrQ2 - sum(pairLoadNewRxy$SrQ2)
              pairLoadNewRxy$WrQ1 = (1+1-(as.numeric(pairLoadNewRxy$PaXy)*2))
              pairLoadNewRxy$WrQ2 = (1+1-(as.numeric(pairLoadNewRxy$PcXy)*2))
              pairLoadNewRxy$WrQ1MinusSum = pairLoadNewRxy$WrQ1 - sum(pairLoadNewRxy$WrQ1)
              pairLoadNewRxy$WrQ2MinusSum = pairLoadNewRxy$WrQ2 - sum(pairLoadNewRxy$WrQ2)
              relatednessR = ((mean(pairLoadNewRxy$SrQ1MinusSum/pairLoadNewRxy$WrQ1MinusSum))+(mean(pairLoadNewRxy$SrQ2MinusSum/pairLoadNewRxy$WrQ2MinusSum)))/2
              
              if(!is.na(relatednessR)) {
                if(relatednessR < 0.0625) { degree = "Unrelated"}
                if(relatednessR >= 0.0625 && relatednessR < 0.1875) { degree = "2nd degree"}
                if(relatednessR >= 0.1875 && relatednessR < 0.3126) { degree = "1st degree"}
                if(relatednessR >= 0.3126) { degree = "Same individual/Twins"}
                
                if(relatednessR == 1) {
                  result$HRC = sprintf("%.4f", round(relatednessR,digits=4))
                  result$counts0 = length(which(pairLoadNewRxy$IacADbcBD %in% 0))
                  result$counts4 = "0"
                  result$Relationship = degree
                } else {
                  result$HRC = sprintf("%.4f", round(relatednessR,digits=4))
                  result$counts0 = length(which(pairLoadNewRxy$IacADbcBD %in% 0))
                  result$counts4 = length(which(pairLoadNewRxy$IacADbcBD %in% 4))
                  result$Relationship = degree
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Clean up temporary files
  unlink(list.files(pattern = temp_prefix))
  
  if(verbose) {
    if(is.na(result$HRC)) {
      cat(paste0("SKIPPING due to no overlapping SNPs\t(",line,"/",nrow(combos),")\n"))
    } else {
      cat(paste0(result$Used_SNPs," SNPs\t", result$HRC," HRC\t(",line,"/",nrow(combos),")\n"))
    }
  }
  
  return(result)
}

# Main processing
if(dyads == TRUE) {
  list_samples = c()
  for(i in list.files(pattern = "\\.ped$")) {
    samp = as.character(strsplit(i,"\\.ped"))
    list_samples = c(list_samples,samp)
  }
  if(length(list_samples) < 2) {
    cat("\n\t # ERROR: Less than 2 '.ped' files found in current folder. A minimum of 2 samples is required\n"); quit()
  } else {
    combos = data.frame(t(combn(list_samples,2)),stringsAsFactors = F)
  }
} else {
  combos = read.table(dyads,header = F, stringsAsFactors = F)
}

sampleList = unique(c(combos[,1],combos[,2]))

if(verbose) {
  cat("######\nMatrix of dyads:\n")
  print(combos)
  cat("\nList of unique samples:\n")
  cat(sampleList)
  cat("\n######\n")
  cat(paste0("\t # Using ", threads, " threads for parallel processing\n"))
}

# Set up parallel processing
cl <- makeCluster(threads)
registerDoParallel(cl)

# Process dyads in parallel
results <- foreach(i=1:nrow(combos), .combine=rbind, .packages=c("data.table")) %dopar% {
  process_dyad(i, combos, freqFile, ignoreThresh, verbose, sampleList)
}

# Stop cluster
stopCluster(cl)

# Export results
write.table(results,"TKGWV2_Results.txt",quote = F,row.names = F,col.names = T,sep="\t")

cat(paste0(" # [",Sys.time(),"] All dyads processed\n"))
cat(paste0(" # [",Sys.time(),"] Results exported to TKGWV2_Results.txt\n"))