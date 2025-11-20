#Obtaining necessary package
# setwd("/panfs/jay/groups/16/saonli/baron063/R")

suppressMessages(library(tidyverse))

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least two arguments: if not, return an error
if (length(args)<=1) {
  stop("A path to the data and the file name for your target population need to be provided.", call.=FALSE)
} else if (length(args)>=2) {
  args[1] -> base_location #Place where the files are located
  args_length = length(args)
  # default output file
  args[2] -> simulated_assoc_name_1 #first .qassoc file of interest (target population)
  if(length(args)==2){
    output_name <- "Output_TLPRS"
  } else {
    args[3] -> simulated_bim_name_1 #first .bim file
    args[4] -> sum_stat_1_name #This is the desired name of the target populations summary stats
    args[5] -> simulated_assoc_name_2 
    args[6] -> simulated_bim_name_2  
    args[7] -> sum_stat_2_name #This is the desired name of the training populations summary stats
  }
}

{
  wd=base_location
  print(wd)
  print(simulated_assoc_name_1)
  print(simulated_bim_name_1)
  print(sum_stat_1_name)
  print(simulated_assoc_name_2)
  print(simulated_bim_name_2)
  print(sum_stat_2_name)
}

#Back to where the data is
setwd(wd) #Argument 1


##testing population data
first_assoc_dat = read.csv(file =simulated_assoc_name_1, 
                        header=TRUE, sep= "") #Argument 2

first_bim_dat = read.csv(file =simulated_bim_name_1, 
                      header = FALSE, sep= "") #Argument 3

first_assoc_data = first_assoc_dat %>% select(SNP, n_eff, BETA, P, A1) 
colnames(first_assoc_data) = c("SNP", "N", "beta", "p", "A1")

first_bim_dat = first_bim_dat %>% select(V5)
colnames(first_bim_dat) = c("A1")

first_sumstats = first_assoc_data

first_sumstats_file = na.omit(first_sumstats %>% select(c(SNP, A1, beta, N, p)))

write.table(first_sumstats_file, file =sum_stat_1_name, 
            sep = " ", quote = FALSE, row.names = FALSE) #Argument 4

print("Wrote testing population summary stats file")
##training data
training_assoc_dat = read.csv(file=simulated_assoc_name_2, 
                        header=TRUE, sep= "") #Argument 5

training_pop_bim_dat = read.csv(file =simulated_bim_name_2,
                      header = FALSE, sep= "") #Argument 6

training_assoc_data = training_assoc_dat %>% select(SNP, n_eff, BETA, P, A1) 
colnames(training_assoc_data) = c("SNP", "N", "beta", "p", "A1")

training_pop_bim_dat = training_pop_bim_dat %>% select(V5)
colnames(training_pop_bim_dat) = c("A1")

training_pop_sumstats = training_assoc_data

training_pop_sumstats_file = na.omit(training_pop_sumstats %>% select(c(SNP, A1, beta, N, p)))

training_pop_sumstats_file_real = training_pop_sumstats_file %>% select(SNP, A1, beta)

write.table(training_pop_sumstats_file_real, file =sum_stat_2_name, 
            sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE) #Argument 7 

print("Wrote training population summary stats file")


