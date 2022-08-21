# This file will generate the plot about ARE between T_{tKW} and T_{KW}. 
# We consider 2 cases here. 
# The result figure is called "ARE_ANOVA.pdf"


alpha = 1:5;

#case1
#Set the parameters
K = 5; theta0 = (1:100)*0.01; R1 = theta0*0;

for(kk in 1:100)
{
theta = numeric(K) + theta0[kk];

#To calculate the deltas
delta = numeric(K - 1);
for(i in 1:(K - 1)){
   for(j in 1:i){
	ins = theta*(alpha/(alpha + alpha[j]) - 1/2);
	insj = theta[j]*sum(ins);
	delta[i] = delta[i] + insj;
    }
    ins = theta*(alpha/(alpha + alpha[i+1]) - 1/2);
    insj = theta[i+1]*sum(ins);
    delta[i] = delta[i] - i*insj;
}

# To calculate the variance for each term
varnum = numeric(K - 1); varden = numeric(K - 1);
thetam = mean(theta);
for(i in 1:(K - 1)){
    varnum[i] = thetam^2*i*(i+1)*(K*(1 - thetam) + K/3);
    varden[i] = i*(i+1)*(K*(1 - thetam) + thetam^2*K/3);
   }
# Calculate the efficiency
R1[kk] = sum(delta^2/varnum)/sum(delta^2/varden)
}


#case2
#Set the parameters
K = 5; thetao = c(0.2, 0.15, 0.3, 0.1, 0.25);
s = (0:70)/100; R2 = s*0;

for(kk in 1:length(s))
{
theta = thetao + s[kk];

#To calculate the deltas
delta = numeric(K - 1);
for(i in 1:(K - 1)){
   for(j in 1:i){
	ins = theta*(alpha/(alpha + alpha[j]) - 1/2);
	insj = theta[j]*sum(ins);
	delta[i] = delta[i] + insj;
    }
    ins = theta*(alpha/(alpha + alpha[i+1]) - 1/2);
    insj = theta[i+1]*sum(ins);
    delta[i] = delta[i] - i*insj;
}

#To calculate the difference of thetas
thetaK = max(theta);
dtheta = numeric(K - 1);
for(i in 1:(K - 1)){
    dtheta[i] = i*theta[i+1] - sum(theta[1:i])/2;
}
dtheta = dtheta*K;

# To calculate the variance for each term
varnum = numeric(K - 1); varden = numeric(K - 1);
thetam = mean(theta);
for(i in 1:(K - 1)){
    varnum[i] = thetam^2*i*(i+1)*(K*(1 - thetam) + K/3);
    varden[i] = i*(i+1)*(K*(1 - thetam) + thetam^2*K/3);
   }

# Calculate the efficiency
R2[kk] = sum((delta + dtheta*thetaK)^2/varnum)/sum((delta + dtheta)^2/varden)
}


### Draw the Plot #####

pdf(file = "ARE_ANOVA.pdf", width = 8, height = 4)
par(mfrow = c(1, 2))
 plot((11:100)/100, R1[11:100], type = "l", lwd = 2, xlab = expression(paste("5 groups, ", theta)), ylab = expression("ARE(T"[tKW]*", T"[KW]*")"))
text(1,58,"(a)", pos = 2, cex = 1.2)
abline(h = 1, lty = 2)
# text(0.15, 3, "y=1", cex = 0.8)

plot((0:70)/100 + 0.3, R2, type = "l", lwd = 2, ylim = c(1, 1.95), xlab = expression(paste("5 groups, ", theta[M])), ylab = expression("ARE(T"[tKW]*", T"[KW]*")"))
text(1,1.9,"(b)", pos = 2, cex = 1.2)
abline(h = 1, lty = 2)
# text(0.32, 1.03, "y=1", cex = 0.8)
dev.off()