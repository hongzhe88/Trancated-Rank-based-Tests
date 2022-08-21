#setwd("D:/dropbox/Dropbox/Penn/NonparaTest")
#setwd("C:/Users/wanjiew/Dropbox/Penn/NonparaTest") # Set the working directory
source('Wilcoxon.zeros.R') # Load the 
require(xlsx)

# read in data set
disease <- read.delim("G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Disease.xls", stringsAsFactors=FALSE)
normal <- read.delim("G_Remove_unclassfied_Renormalized_Merge_Rel_MetaPhlAn_Result_Normal.xls", stringsAsFactors=FALSE)
rownames(disease) <- disease[,"X"]
rownames(normal) <- normal[,"X"]

# select observations with the same name
intersection <- intersect(normal[,"X"], disease[,"X"])
disease <- disease[intersection,]
normal <- normal[intersection,]

# select the disease data at time point 01
disease <- disease[, grep(".01", colnames(disease), value = T,fixed = T)]
# get rid of the "X" column
normal[,"X"] = NULL
# Normalize the genera again
disease= t(t(disease)/colSums(disease))*100;
normal= t(t(normal)/colSums(normal))*100;
ddata <- as.matrix(disease)
ndata <- as.matrix(normal)

# remove data containing only zeros
ddata <- ddata[(rowSums(disease) != 0) & (rowSums(normal) != 0),]
ndata <- ndata[(rowSums(disease) != 0 & rowSums(normal) != 0),]

num <- dim(ddata)[1] #number of observations
pvalues <- matrix(NA, nrow = num, ncol = 4)
colnames(pvalues) <- c("newStat_perm", "newStat_chisq", "Wilcoxon", "twopart")  
rownames(pvalues) <- rownames(ddata)


#G1 is disease while G2 is normal
perc1 = perc2 = numeric(num);
for (i in (1:num)){
  G1 <- ddata[i,]
  G2 <- ndata[i,]
  perc1[i] <- mean(G1==0); n1 = length(G1)
  perc2[i] <- mean(G2==0); n2 = length(G2)
  pvalues[i, "newStat_perm"] <- Wilcoxon.zeros(G1, G2, perm = TRUE)$p.value         #pvalue for the new statistic
#  pvalues[i, "newStat_chisq"] <- Wilcoxon.zeros(G1, G2)$p.value       #pvalue for the new statistic
  pvalues[i, "Wilcoxon"] <- wilcox.test(G1, G2, paired = F, exact = FALSE)$p.value #pvalue for the two sample Ttest ("two.sided")
  G1score = G1[G1!=0]; G2score = G2[G2!=0]; 
  stat1 = (perc1[i] - perc2[i])^2/(perc1[i]*(1-perc1[i])/n1+perc2[i]*(1-perc2[i])/n2);
  n1 = length(G1score); n2 = length(G2score); r = rank(c(G1score, G2score), tie ="average");
  stat2 = (sum(r[1:n1]) - n1*(n1+n2+1)/2)/sqrt(n1*n2*(n1+n2+1)/12)
  pvalues[i, "twopart"] <- pchisq(stat1+stat2^2, 2, lower.tail = FALSE) #pvalue for the two sample Ttest ("two.sided")
}
perc1 = perc1[!is.na(rowSums(pvalues))]; perc2 = perc2[!is.na(rowSums(pvalues))];
pvalues = pvalues[!is.na(rowSums(pvalues)),]
View(pvalues)

# filename <- file.path(getwd(), paste0("genus__disease_different_time_pvalues.xlsx"))
# write.xlsx(pvalues,filename)
alpha <- 0.05
pwr <- colSums(pvalues < alpha)/num
sigGenius_Perm <- rownames(pvalues[pvalues[,"newStat_perm"] < alpha,])
length(sigGenius_Perm ) #26
sigGenius_Chisq <- rownames(pvalues[pvalues[,"newStat_chisq"] < alpha,])
length(sigGenius_Chisq )  #23
sigGenius_Wilcox <- rownames(pvalues[pvalues[,"Wilcoxon"] < alpha,])
length(sigGenius_Wilcox )   #20
sigGenius_twopart <- rownames(pvalues[pvalues[,"twopart"] < alpha,])
length(sigGenius_twopart )   #30

dif = c(setdiff(sigGenius_Chisq, sigGenius_Wilcox), setdiff(sigGenius_Wilcox, sigGenius_Chisq));
difvalue = cbind(ndata[dif,], ddata[dif,])
nn = dim(ndata)[2]; nd = dim(ddata)[2]
label = c(rep("Normal", nn), rep("Crohn's Disease", nd))
colnames(difvalue) <- label;

dif2 = c(setdiff(sigGenius_Chisq, sigGenius_twopart), setdiff(sigGenius_twopart, sigGenius_Chisq));
difvalue2 = cbind(ndata[dif2,], ddata[dif2,])
colnames(difvalue2) <- label;


pdf(file = "Wilcoxon_5plot.pdf", height = 7, width = 12)
par(mfrow = c(2,3), mai = c(.7, 1,.7,.2))
for(i in 1:5){
  main1 = paste(strsplit(dif[i], "_")[[1]][3], ", ", sep = "");
  main2 = paste("=", substr(as.character(round(pvalues[dif[i], "newStat_chisq"], 3)), 2, 5), ", ", sep = "");
  main3 = paste("=", substr(as.character(round(pvalues[dif[i], "Wilcoxon"], 3)), 2, 5), sep = "");
  x1 = difvalue[i, 1:26]; x2 = difvalue[i, 27:111];
  perc1 = mean(x1 == 0); perc2 = mean(x2 == 0);
  x1 = log(x1[x1 != 0]); x2 = log(x2[x2 != 0]); 
  label = c(rep(paste("Normal (", round(perc1*100), "%)", sep = ""), length(x1)), 
            rep(paste("Crohn (", round(perc2*100), "%)", sep = ""), length(x2)))
  if(i < 5) {boxplot(c(x1, x2)~label, col = "white", cex.axis = 1.7, xlab = NULL, ylab = NULL); ymax = max(x1, x2);}
  else {boxplot(c(x1, x2)~label, col = "white", 
                cex.axis = 1.7, xlab = NULL, ylab = NULL); ymax = 8.5;}
  title(main = bquote(.(main1) ~ p[s] ~ .(main2) ~ p[S] ~ .(main3)), cex.main=1.8);
  if(i == 1 | i == 4 |i==7){
    title(ylab = "Relative Abundance (log scale)", cex.lab = 2)
  }
}
dev.off()

#save.image(file = "Normal_vs_Disease.Rdata")
