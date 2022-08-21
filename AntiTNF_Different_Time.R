# setwd("C:/Users/wanjiew/Box Sync/NonparaTest")
source('ANOVA.zeros.R')
require(xlsx)

# read in data set
disease <- read.delim("G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Disease.xls", stringsAsFactors=FALSE)
clinic  <- read.csv("COMBO_PLEASE_Sample_Information_Taxa.csv")
rownames(disease) <- disease[, "X"]
disease[,"X"] = NULL;

# Find the samples with disease in clinic
clinic <- clinic[!(is.na(clinic[,"Treatment"])), ]

# Find the overlapped samples in clinic and disease
clinic_obs <- clinic[, "X"];
clinic_obs1 <- substr(clinic_obs, start = 1, stop = 4)
clinic_obs2 <- substr(clinic_obs, start = 6, stop = 7)
clinic_obs <- paste("X", clinic_obs1, ".", clinic_obs2, sep = "");
intersection <- intersect(clinic_obs, colnames(disease))
disease <- disease[,intersection]

# select antiTNF samples
indtnf = which(clinic[,"Treatment"] == "antiTNF");
disease = disease[ , indtnf+1];
clinic = clinic[indtnf, ]

# Table to save p-values
num <- dim(disease)[1] #number of generals
pvalues <- matrix(NA, nrow = num, ncol = 3)
colnames(pvalues) <- c("newStat_perm", "newStat_chisq", "Friedman")  
rownames(pvalues) <- rownames(disease)

# Generate the table at different time points
sample = dim(disease)[2]; Nsample = numeric(num);
clinic[, "Subject"] <- as.factor(clinic[, "Subject"])
for (i in (1:num)){

  data = matrix(NA, ncol = 4, nrow = nlevels(clinic[, "Subject"]));
  rownames(data) <- levels(clinic[, "Subject"])
  for (j in 1:sample){
    data[clinic[j, "Subject"], clinic[j, "Time"]] = disease[i, j];
  }
  data <- data[(rowSums(data) != 0) & (!is.na(rowSums(data))), ]

  Nsample[i] = length(data)/4;
  if(length(data) <= 4) next;
  data.test = list(data[,1], data[,2], data[,3], data[,4]);

  w = ANOVA.zeros(data.test)$statistics; permw = numeric(10000);
  for(ii in 1:10000){
	data_perm = data;
	for(jj in 1:nrow(data)) data_perm[jj,] = sample(data[jj,], 4);
	data_perm = list(data_perm[,1], data_perm[,2], data_perm[,3], data_perm[,4]);
	permw[ii] = ANOVA.zeros(data_perm)$statistics
  }
  permw = as.numeric(permw);
  pvalues[i, "newStat_perm"] <- sum(abs(w) < abs(permw))/10000;         #pvalue for the new statistic
  pvalues[i, "newStat_chisq"] <- ANOVA.zeros(data.test)$p.value       #pvalue for the new statistic
  pvalues[i, "Friedman"] <- friedman.test(data)$p.value #pvalue for Repeated Measurement Tests
  print(i)
}
pvalues = pvalues[!is.na(rowSums(pvalues)),];
View(pvalues)

#filename <- file.path(getwd(), paste0("genus_normal_vs_disease_pvalues.xlsx"))
#write.xlsx(pvalues,filename)

alpha <- 0.05
pwr <- colSums(pvalues < alpha)/num
sigGenius_Perm <- rownames(pvalues[pvalues[,"newStat_perm"] < alpha,])
length(sigGenius_Perm ) #7
sigGenius_Chisq <- rownames(pvalues[pvalues[,"newStat_chisq"] < alpha,])
length(sigGenius_Chisq )  #2
sigGenius_Friedman<- rownames(pvalues[pvalues[,"Friedman"] < alpha,])
length(sigGenius_Friedman)   #3



setdiff(sigGenius_Friedman, sigGenius_Chisq)
setdiff(sigGenius_Perm, sigGenius_Friedman)
setdiff(sigGenius_Perm, sigGenius_Chisq)

dif = c(setdiff(sigGenius_Perm, sigGenius_Friedman), setdiff(sigGenius_Friedman, sigGenius_Perm));
difvalue = disease[dif,];
label = rep(c("Baseline", "Week 1", "Week 4", "Week 8"), ncol(disease)/4)
colnames(difvalue) <- label;

pdf(file = "RMtest_new.pdf", width = 11, height = 6)
par(mfcol = c(2, 4), mai = c(.55, 1, .35,.2))
for(i in 1:4){
  data = matrix(NA, ncol = 4, nrow = nlevels(clinic[, "Subject"]));
  rownames(data) <- levels(clinic[, "Subject"])
  colnames(data) <- c("0", "1", "4", "8")
  for (j in 1:sample)
    data[clinic[j, "Subject"], clinic[j, "Time"]] = disease[dif[i], j];
  data = data[(rowSums(data)!=0) & !(is.na(rowSums(data))), ];
  main1 = paste(strsplit(dif[i], "_")[[1]][3], ",", sep = "");
  main2 = paste("=", substr(round(pvalues[dif[i], "newStat_perm"], 3), 2, 5), ",", sep = "");
  main3 = paste("=", substr(round(pvalues[dif[i], "Friedman"], 3), 2, 5),sep = "");
  data2 = apply(data == 0, 2, mean)
#  barplot(data2, main = bquote(.(main1) ~ p[kw] ~ .(main2) ~ p[KW] ~ .(main3)), cex.names = 2, cex.main = 1.8)
  barplot(data2, xlab = strsplit(dif[i], "_")[[1]][3], cex.names = 2, cex.axis = 1.8, cex.lab = 1.8)
  if(i == 1)
    title(ylab = "Percentage of Zeroes", cex.lab = 2)
  
  nz = apply(data != 0, 2, sum);
  label = c(rep("0", nz[1]), rep("1", nz[2]), rep("4", nz[3]), rep("8", nz[4])); 
  data1 = c(data[data[,1]!=0, 1], data[data[,2]!=0, 2], data[data[,3]!=0, 3], data[data[,4]!=0, 4]);
#  boxplot(data1~label, main = bquote(.(main1) ~ p[kw] ~ .(main2) ~ p[KW] ~ .(main3)), cex.main = 1.8, cex.axis = 2)
  boxplot(log(data1)~label, main = bquote(p[kw] ~ .(main2) ~ p[KW] ~ .(main3)), cex.main = 1.8, cex.axis = 2, 
          ylab = NULL, xlab = NULL)
  if(i == 1)
    title(ylab = "Relative Abundance (log scale)", cex.lab = 2)
  print(nrow(data))
}
dev.off()