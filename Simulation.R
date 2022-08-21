## This file is to run the simulation experiments ###
### It will generate 7 Rdata files, where the first file is for Experiment 1, 
## the 2-7 files are for Experiment 2
## It will generate "power.pdf", which corresponds to Figure 4 in the paper

#############################################################################
#Experiment 1. Size under null hypothesis##########################################
#############################################################################

rep = 100000;

theta = 0.5;
alpha = 2; beta = 2;

Narray = c(30, 60, 100, 300, 600, 900);
pkw = array(0, dim = c(3, length(Narray), rep));
pw = array(0, dim = c(3, length(Narray), rep));
sizekw = array(0, dim = c(3, 3, length(Narray)));
sizew = array(0, dim = c(3, 3, length(Narray)));
w = pw;

jj = 1; K = 2; Nprop = c(.65, 1);
for(ii in 1:length(Narray)){
   Nbase = Narray[ii];
   for(kk in 1:rep){
	N = round(Nbase*Nprop);
	###########Generate Data ###################
	Q = rbinom(sum(N), 1, theta);
	X1 = rbeta(N[1], alpha, beta)
	X1 = Q[1:N[1]]*X1;
	X2 = rbeta(N[2], alpha, beta)
	X2 = Q[N[1] + 1:N[2]]*X2;

	###########Wilcoxon Rank Sum Test ############
	pkw[jj, ii, kk] = wilcox.test(X1, X2)$p.value;
	
	###########Our Test ###################
	n = c(sum(Q[1:N[1]]), sum(Q[N[1] + 1:N[2]]))
	prop = n/N; pmax = max(prop); pmean = mean(prop);
	Ntrun = round(pmax*N);
	X1trun = c(X1[Q[1:N[1]] == 1], rep(0, Ntrun[1] - n[1]));
	X2trun = c(X2[Q[N[1] + 1:N[2]] == 1], rep(0, Ntrun[2]- n[2]));
	rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun));
	r = sum(rankdata[1:Ntrun[1]]);
	s = r - Ntrun[1]*(sum(Ntrun) + 1)/2;

	v1 = pmean*(1 - pmean);
	var1 = N[1]*N[2]*pmean*v1*(sum(N)*pmean + 3*(1 - pmean)/2) + sum(N^2)*v1^2*5/4 + 2*sqrt(2/pi)*pmean*(sum(N)*v1)^(1.5)*sqrt(N[1]*N[2]);
	var1 = var1/4;
	var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
	vars = var1 + var2;
	w[jj, ii, kk] = s/sqrt(vars);
	}
	pw[jj,ii,] = 2*(1 - pnorm(abs(w[jj,ii,])));
	print(ii)
}

for(ii in 1:length(Narray)){
	sizekw[jj, 1, ii] = mean(pkw[jj, ii, ] < 0.05);
	sizew[jj, 1, ii] = mean(pw[jj, ii, ] < 0.05);
	sizekw[jj, 2, ii] = mean(pkw[jj, ii, ] < 0.01);
	sizew[jj, 2, ii] = mean(pw[jj, ii, ] < 0.01);
	sizekw[jj, 3, ii] = mean(pkw[jj, ii, ] < 0.001);
	sizew[jj, 3, ii] = mean(pw[jj, ii, ] < 0.001);
}

jj = 2; K = 3;
for(ii in 1:length(Narray)){
   N  = Narray[ii];
   for(kk in 1:rep){

	###########Generate Data ###################
	Q = rbinom(3*N, 1, theta);
	X = rbeta(3*N, alpha, beta)
	X = Q*X;
	X1 = X[1:N];
	X2 = X[N + 1:N];
	X3 = X[N*2 + 1:N];

	###########Kruskal Wallis Test ############
	data = list(X1, X2, X3);
	pkw[jj, ii, kk] = kruskal.test(data)$p.value;
	
	###########Our Test ###################
	n = c(sum(Q[1:N]), sum(Q[N+1:N]), sum(Q[N*2+1:N])); nmax = max(n);
	ind = which(Q == 1);
	X1trun = c(X[ind[ind <= N]], rep(0, nmax - n[1]));
	X2trun = c(X[ind[ind > N & ind <= 2*N]], rep(0, nmax - n[2]));
	X3trun = c(X[ind[ind > 2*N]], rep(0, nmax - n[3]));
	rankdata = K*nmax + 1 - rank(c(X1trun, X2trun, X3trun));
	r = c(sum(rankdata[1:nmax]), sum(rankdata[nmax+1:nmax]), sum(rankdata[2*nmax+1:nmax]) ); 
	s = r - nmax*(K*nmax + 1)/2;
	u = c(s[1] - s[2], sum(s[1:2]) - 2*s[3]);
	thetam = mean(n)/N;
	
	simun = matrix(rbinom(5000*K, N, thetam), nrow = 5000, ncol = K);
	simunmax = apply(simun, 1, max); simusum = apply(simun, 1, sum);
	varsimu1 = numeric(K - 1); varsimu2 = varsimu1;
	varsimu1[1] = mean(simunmax^2*(simun[,1] - simun[,2])^2)*K^2/4;
	varsimu2[1] = mean((simun[,1]*simun[,2]*4 + simun[,3]*(simun[,1] + simun[,2]))*(simusum+1)/12);
	for(ss in 2:(K-1)){
		simuss = apply(simun[,1:ss], 1, sum);
		varsimu1[ss] = mean(simunmax^2*(simuss - simun[,ss+1]*ss)^2)*K^2/4;
		varsimu2[ss] = mean((simuss*simun[, ss+1]*(ss+1)^2 + (simusum - simuss - simun[,ss+1])*(simuss + ss^2*simun[,ss+1]))*(simusum+1)/12);
	}
	varsimu = varsimu1 + varsimu2;
	w[jj, ii, kk] = sum(u^2/varsimu);
	}
	pw[jj,ii,] = pchisq(w[jj,ii,], 2, lower.tail = F);
	print(ii)
}

for(ii in 1:length(Narray)){
	sizekw[jj, 1, ii] = mean(pkw[jj, ii, ] < 0.05);
	sizew[jj, 1, ii] = mean(pw[jj, ii, ] < 0.05);
	sizekw[jj, 2, ii] = mean(pkw[jj, ii, ] < 0.01);
	sizew[jj, 2, ii] = mean(pw[jj, ii, ] < 0.01);
	sizekw[jj, 3, ii] = mean(pkw[jj, ii, ] < 0.001);
	sizew[jj, 3, ii] = mean(pw[jj, ii, ] < 0.001);
}


jj = 3; K = 3; Nprop = c(0.7, 1, 1.5);
for(ii in 1:length(Narray)){
   Nbase = Narray[ii];
   for(kk in 1:rep){
	N = round(Nbase*Nprop);
	###########Generate Data ###################
	Q = rbinom(sum(N), 1, theta);
	X1 = rbeta(N[1], alpha, beta)
	X1 = Q[1:N[1]]*X1;
	X2 = rbeta(N[2], alpha, beta)
	X2 = Q[N[1] + 1:N[2]]*X2;
	X3 = rbeta(N[3], alpha, beta)
	X3 = Q[N[1] + N[2] +1:N[3]]*X3;

	###########Kruskal Wallis Test ############
	data = list(X1, X2, X3);
	pkw[jj, ii, kk] = kruskal.test(data)$p.value;
	
	###########Our Test ###################
	n = c(sum(Q[1:N[1]]), sum(Q[N[1] + 1:N[2]]), sum(Q[N[1] + N[2] +1:N[3]]))
	prop = n/N; pmax = max(prop);
	Ntrun = round(pmax*N);
	X1trun = c(X1[Q[1:N[1]] == 1], rep(0, Ntrun[1] - n[1]));
	X2trun = c(X2[Q[N[1] + 1:N[2]] == 1], rep(0, Ntrun[2]- n[2]));
	X3trun = c(X3[Q[N[1] + N[2] +1:N[3]] == 1], rep(0, Ntrun[3]- n[3]));
	rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun, X3trun));
	r = c(sum(rankdata[1:Ntrun[1]]), sum(rankdata[Ntrun[1]+1:Ntrun[2]]), sum(rankdata[sum(Ntrun[1:2])+1:Ntrun[3]]) ); 
	s = r - Ntrun*(sum(Ntrun) + 1)/2;
	u = c(N[2]*s[1] - N[1]*s[2], N[3]*sum(s[1:2]) - s[3]*sum(N[1:2]));
	u = u/Nbase^2;

	thetam = mean(prop);
	simun = matrix(0, nrow = 5000, ncol = K); simup = simun;
	for(ss in 1:K){
		simun[,ss] = rbinom(5000, N[ss], thetam);		
		simup[,ss] = simun[,ss]/N[ss];
	}
	simupmax = apply(simup, 1, max);
	varsimu = numeric(K - 1);
	varsimu[1] = Nprop[2]^2*mean(simupmax^2*(simup[,1] - simup[,2])^2)*N[1]^2;
	for(ss in 2:(K-1)){
		varsimu[ss] = Nprop[ss+1]^2*mean(simupmax^2*(apply(simun[,1:ss], 1, sum) - simup[,ss+1]*sum(N[1:ss]))^2);
	}
	varsimu = varsimu*(sum(Nprop))^2/4;

	varu2 = c(Nprop[1]*Nprop[2], Nprop[3]*sum(Nprop))*(Nprop[1] + Nprop[2])*sum(Nprop)*thetam^2*(sum(N)*thetam + 1)/12;
	varu = varsimu + varu2;
	w[jj, ii, kk] = sum(u^2/varu);
	}
	pw[jj,ii,] = pchisq(w[jj,ii,], 2, lower.tail = F);
	print(ii)
}

for(ii in 1:length(Narray)){
	sizekw[jj, 1, ii] = mean(pkw[jj, ii, ] < 0.05);
	sizew[jj, 1, ii] = mean(pw[jj, ii, ] < 0.05);
	sizekw[jj, 2, ii] = mean(pkw[jj, ii, ] < 0.01);
	sizew[jj, 2, ii] = mean(pw[jj, ii, ] < 0.01);
	sizekw[jj, 3, ii] = mean(pkw[jj, ii, ] < 0.001);
	sizew[jj, 3, ii] = mean(pw[jj, ii, ] < 0.001);
}
save(w, pkw, pw, sizew, sizekw, Narray, theta, alpha, beta, rep, file = "ANOVAexp1.RData")
rm(list = ls())


#############################################################################
#Experiment 2a-b: equal size tests with the same theta and different nonzero functions
#############################################################################

rep = 10000;
K = 3;

theta = c(0.5, 0.5, 0.5);
alphamatrix = matrix(c(1.5, 2, 2.5, 1.8, 2, 2.2, 1.5, 2, 2.5, 1.8, 2, 2.2), nrow = 4, ncol = K, byrow = T);
betamatrix = matrix(c(2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4), nrow = 4, ncol = K, byrow = T);
pkw = array(0, dim = c(2, 20, rep));
pw = array(0, dim = c(2, 20, rep));
powerkw = matrix(0, nrow = 2, ncol = 20);
powerw = matrix(0, nrow = 2, ncol = 20);
w = pw;

for(ii in 1:20){
   N = ii*100;
   for(jj in 1:2){
      alpha = alphamatrix[jj,];
    	beta = betamatrix[jj, ];
   for(kk in 1:rep){
     set.seed(20220*ii+10000*jj+kk)
     ###########Generate Data ###################
     Q = rbinom(3*N, 1, mean(theta));
     X1 = rbeta(N, alpha[1], beta[1])
     X1 = Q[1:N]*X1;
     X2 = rbeta(N, alpha[2], beta[2])
     X2 = Q[(N + 1):(2*N)]*X2;
     X3 = rbeta(N, alpha[3], beta[3])
     X3 = Q[(2*N + 1):(3*N)]*X3;
     
     
     ###########Kruskal Wallis Test ############
     data = list(X1, X2, X3);
     pkw[jj, ii, kk] = kruskal.test(data)$p.value;
     
     ###########Our Test ###################
     n = c(sum(Q[1:N]), sum(Q[(N + 1):(2*N)]), sum(Q[(2*N + 1):(3*N)])); nmax = max(n);
     X1trun = c(X1[Q[1:N] == 1], rep(0, nmax - n[1]));
     X2trun = c(X2[Q[(N + 1):(2*N)] == 1], rep(0, nmax - n[2]));
     X3trun = c(X3[Q[(2*N + 1):(3*N)] == 1], rep(0, nmax - n[3]));
     rankdata = K*nmax + 1 - rank(c(X1trun, X2trun, X3trun));
     r = c(sum(rankdata[1:nmax]), sum(rankdata[nmax+1:nmax]), sum(rankdata[2*nmax+1:nmax]) );
     s = r - nmax*(K*nmax + 1)/2;
     u = c(s[1] - s[2], sum(s[1:2]) - 2*s[3])/N^2;
     
     thetam = mean(n)/N;
     simun = matrix(rbinom(5000*K, N, thetam), nrow = 5000, ncol = K);
     simunmax = apply(simun, 1, max); simusum = apply(simun, 1, sum);
     varsimu1 = numeric(K - 1); varsimu2 = varsimu1;
     varsimu1[1] = mean(simunmax^2*(simun[,1] - simun[,2])^2)*K^2/4;
     varsimu2[1] = mean((simun[,1]*simun[,2]*4 + simun[,3]*(simun[,1] + simun[,2]))*(simusum+1)/12);
     for(ss in 2:(K-1)){
       simuss = apply(simun[,1:ss], 1, sum);
       varsimu1[ss] = mean(simunmax^2*(simuss - simun[,ss+1]*ss)^2)*K^2/4;
       varsimu2[ss] = mean((simuss*simun[, ss+1]*(ss+1)^2 + (simusum - simuss - simun[,ss+1])*(simuss + ss^2*simun[,ss+1]))*(simusum+1)/12);
     }
     varsimu = varsimu1 + varsimu2;
     w[jj, ii, kk] = sum(u^2/varsimu);
   }
	pw = pchisq(w, 2, lower.tail = F);
	powerkw[jj, ii] = mean(pkw[jj, ii, ] < 0.001);
	powerw[jj, ii] = mean(pw[jj, ii, ] < 0.001);
   }
   print(ii)
}

save(pw, pkw, powerkw, powerw, alphamatrix,N, betamatrix, theta, K, rep, file = "ANOVAexp2.RData")


#################################################################################
#Experiment 2c-d: equal size tests with the same nonzero functions and different theta#
##############################################################################
rep = 10000;
K = 3;

thetamatrix = c(0.1, 0.15, 0.22);
alpha = 1/2; beta = 2;
pkw = array(0, dim = c(2, 15, rep));
pw = array(0, dim = c(2, 15, rep));
powerkw = matrix(0, nrow = 2, ncol = 15);
powerw = matrix(0, nrow = 2, ncol = 15);
kw = pkw; w = pw;

for(ii in 1:15){
  N = ii*20;
  for(jj in 1:2){
    theta = thetamatrix + 0.05*(jj-1);
    for(kk in 1:rep){
      set.seed(20220*ii+10000*jj+kk)
      ###########Generate Data ###################
      X = rbeta(3*N, alpha, beta);
      Q1 = rbinom(N, 1, theta[1]);
      X1 = X[1:N]*Q1;
      Q2 = rbinom(N, 1, theta[2]);
      X2 = X[(N + 1):(2*N)]*Q2;
      Q3 = rbinom(N, 1, theta[3]);
      X3 = X[(2*N + 1):(3*N)]*Q3;
      
      ###########Kruskal Wallis Test ############
      data = list(X1, X2, X3);
      pkw[jj, ii, kk] = kruskal.test(data)$p.value;
      
      ###########Our Test ###################
      n = c(sum(Q1), sum(Q2), sum(Q3)); nmax = max(n);
      X1trun = c(X1[Q1 == 1], rep(0, nmax - n[1]));
      X2trun = c(X2[Q2 == 1], rep(0, nmax - n[2]));
      X3trun = c(X3[Q3 == 1], rep(0, nmax - n[3]));
      rankdata = K*nmax + 1 - rank(c(X1trun, X2trun, X3trun));
      r = c(sum(rankdata[1:nmax]), sum(rankdata[nmax+1:nmax]), sum(rankdata[2*nmax+1:nmax]) ); 
      s = r - nmax*(K*nmax + 1)/2;
      u = c(s[1] - s[2], sum(s[1:2]) - 2*s[3]);
      thetam = mean(n)/N;
      
      simun = matrix(rbinom(5000*K, N, thetam), nrow = 5000, ncol = K);
      simunmax = apply(simun, 1, max); simusum = apply(simun, 1, sum);
      varsimu1 = numeric(K - 1); varsimu2 = varsimu1;
      varsimu1[1] = mean(simunmax^2*(simun[,1] - simun[,2])^2)*K^2/4;
      varsimu2[1] = mean((simun[,1]*simun[,2]*4 + simun[,3]*(simun[,1] + simun[,2]))*(simusum+1)/12);
      for(ss in 2:(K-1)){
        simuss = apply(simun[,1:ss], 1, sum);
        varsimu1[ss] = mean(simunmax^2*(simuss - simun[,ss+1]*ss)^2)*K^2/4;
        varsimu2[ss] = mean((simuss*simun[, ss+1]*(ss+1)^2 + (simusum - simuss - simun[,ss+1])*(simuss + ss^2*simun[,ss+1]))*(simusum+1)/12);
      }
      varsimu = varsimu1 + varsimu2;
      w[jj, ii, kk] = sum(u^2/varsimu);
    }
    pw = pchisq(w, 2, lower.tail = F);
    powerkw[jj, ii] = mean(pkw[jj, ii, ] < 0.001);
    powerw[jj, ii] = mean(pw[jj, ii, ] < 0.001);
  }
  print(ii)
}

save(w, pkw, pw, powerkw, powerw, alpha, beta,N, thetamatrix, rep, file = "ANOVAexp3.RData")
rm(list = ls())


#############################################################################
#Experiment 2e-f, unequal size tests with the different nonzero functions and different theta
#############################################################################

rep = 10000;
K = 3;

theta = c(0.4, 0.5, 0.6); Nprop = c(0.8, 1, 1.5);
alphamatrix = matrix(c(1.5, 2, 2.5, 1.7, 2, 2.3, 1.5, 2, 2.5, 1.7, 2, 2.3), nrow = 4, ncol = K, byrow = T);
betamatrix = matrix(c(2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4), nrow = 4, ncol = K, byrow = T);
pkw = array(0, dim = c(2, 20, rep));
pw = array(0, dim = c(2, 20, rep));
powerkw = matrix(0, nrow = 2, ncol = 20);
powerw = matrix(0, nrow = 2, ncol = 20);
w = pw;

for(ii in 1:20){
  Nbase = ii*100;
  for(jj in 1:2){
    alpha = alphamatrix[jj,];
    beta = betamatrix[jj, ];
    for(kk in 1:rep){
      set.seed(20220*ii+10000*jj+kk)
      N = round(Nbase*Nprop);
      ###########Generate Data ###################
      Q = rbinom(sum(N), 1, mean(theta));
      X1 = rbeta(N[1], alpha[1], beta[1])
      X1 = Q[1:N[1]]*X1;
      X2 = rbeta(N[2], alpha[2], beta[2])
      X2 = Q[N[1] + 1:N[2]]*X2;
      X3 = rbeta(N[3], alpha[3], beta[3])
      X3 = Q[N[1] + N[2] +1:N[3]]*X3;
      
      ###########Kruskal Wallis Test ############
      data = list(X1, X2, X3);
      pkw[jj, ii, kk] = kruskal.test(data)$p.value;
      
      ###########Our Test ###################
      n = c(sum(Q[1:N[1]]), sum(Q[N[1] + 1:N[2]]), sum(Q[N[1] + N[2] +1:N[3]]))
      prop = n/N; pmax = max(prop);
      Ntrun = round(pmax*N);
      X1trun = c(X1[Q[1:N[1]] == 1], rep(0, Ntrun[1] - n[1]));
      X2trun = c(X2[Q[N[1] + 1:N[2]] == 1], rep(0, Ntrun[2]- n[2]));
      X3trun = c(X3[Q[N[1] + N[2] +1:N[3]] == 1], rep(0, Ntrun[3]- n[3]));
      rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun, X3trun));
      r = c(sum(rankdata[1:Ntrun[1]]), sum(rankdata[Ntrun[1]+1:Ntrun[2]]), sum(rankdata[sum(Ntrun[1:2])+1:Ntrun[3]]) ); 
      s = r - Ntrun*(sum(Ntrun) + 1)/2;
      u = c(N[2]*s[1] - N[1]*s[2], N[3]*sum(s[1:2]) - s[3]*sum(N[1:2]))/Nbase^2;
      
      thetam = mean(prop);
      simun = matrix(0, nrow = 5000, ncol = K); simup = simun;
      for(ss in 1:K){
        simun[,ss] = rbinom(5000, N[ss], thetam);		
        simup[,ss] = simun[,ss]/N[ss];
      }
      simupmax = apply(simup, 1, max);
      varsimu = numeric(K - 1);
      varsimu[1] = Nprop[2]^2*mean(simupmax^2*(simup[,1] - simup[,2])^2)*N[1]^2;
      for(ss in 2:(K-1)){
        varsimu[ss] = Nprop[ss+1]^2*mean(simupmax^2*(apply(simun[,1:ss], 1, sum) - simup[,ss+1]*sum(N[1:ss]))^2);
      }
      varsimu = varsimu*(sum(Nprop))^2/4;
      
      varu2 = c(Nprop[1]*Nprop[2], Nprop[3]*sum(Nprop))*(Nprop[1] + Nprop[2])*sum(Nprop)*thetam^2*(sum(N)*thetam + 1)/12;
      varu = varsimu + varu2;
      w[jj, ii, kk] = sum(u^2/varu);
    }
    pw = pchisq(w, 2, lower.tail = F);
    powerkw[jj, ii] = mean(pkw[jj, ii, ] < 0.001);
    powerw[jj, ii] = mean(pw[jj, ii, ] < 0.001);
  }
  print(ii)
}

save(pw, pkw, powerkw, powerw, alphamatrix,N, Nprop, betamatrix, theta, K, rep, file = "ANOVAexp4.RData")
rm(list = ls())


###### Set up Sequencing Depth ######
### This file is to check what will happen when the sequencing depth is not enough ###
### We re-run simulation 2, by randomly setting some reads to be 0 ###

#############################################################################
#Exp 2a-b, equal size tests with different nonzero functions and the same theta
#############################################################################

rep = 10000;
K = 3;

theta = c(0.5, 0.5, 0.5);
alphamatrix = matrix(c(1.5, 2, 2.5, 1.8, 2, 2.2, 1.5, 2, 2.5, 1.8, 2, 2.2), nrow = 4, ncol = K, byrow = T);
betamatrix = matrix(c(2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4), nrow = 4, ncol = K, byrow = T);
pkw = array(0, dim = c(2, 20, rep));
pw = array(0, dim = c(2, 20, rep));
powerkw = matrix(0, nrow = 2, ncol = 20);
powerw = matrix(0, nrow = 2, ncol = 20);
w = pw;

for(ii in 1:20){
  N = ii*100;
  for(jj in 1:2){
    alpha = alphamatrix[jj,];
    beta = betamatrix[jj, ];
    for(kk in 1:rep){
      set.seed(20220*ii+10000*jj+kk)
      ###########Generate Data ###################
      Q = rbinom(3*N, 1, mean(theta));
      X1 = rbeta(N, alpha[1], beta[1])
      X1 = Q[1:N]*X1;
      X2 = rbeta(N, alpha[2], beta[2])
      X2 = Q[(N + 1):(2*N)]*X2;
      X3 = rbeta(N, alpha[3], beta[3])
      X3 = Q[(2*N + 1):(3*N)]*X3;
      
      ###some reads are missing ####
      ##when the reads are small, there is 1/2 probability missing ##
      ## small means the value is smaller than a ###
      a = 0.1;
      M = rbinom(length(X1), 1, 0.5);
      M[X1 > a] = 1;
      X1 = X1*M; 
      Q1 = 1*(X1 != 0);
      
      M = rbinom(length(X2), 1, 0.5);
      M[X2 > a] = 1;
      X2 = X2*M; 
      Q2 = 1*(X2 != 0);
      
      M = rbinom(length(X3), 1, 0.5);
      M[X3 > a] = 1;
      X3 = X3*M; 
      Q3 = 1*(X3 != 0);
      
      Q = c(Q1, Q2, Q3);
      ###########Kruskal Wallis Test ############
      data = list(X1, X2, X3);
      pkw[jj, ii, kk] = kruskal.test(data)$p.value;
      
      ###########Our Test ###################
      n = c(sum(Q1), sum(Q2), sum(Q3)); nmax = max(n);
      X1trun = c(X1[Q1 == 1], rep(0, nmax - n[1]));
      X2trun = c(X2[Q2 == 1], rep(0, nmax - n[2]));
      X3trun = c(X3[Q3 == 1], rep(0, nmax - n[3]));
      rankdata = K*nmax + 1 - rank(c(X1trun, X2trun, X3trun));
      r = c(sum(rankdata[1:nmax]), sum(rankdata[nmax+1:nmax]), sum(rankdata[2*nmax+1:nmax]) );
      s = r - nmax*(K*nmax + 1)/2;
      u = c(s[1] - s[2], sum(s[1:2]) - 2*s[3]);
      
      thetam = mean(n)/N;
      simun = matrix(rbinom(5000*K, N, thetam), nrow = 5000, ncol = K);
      simunmax = apply(simun, 1, max); simusum = apply(simun, 1, sum);
      varsimu1 = numeric(K - 1); varsimu2 = varsimu1;
      varsimu1[1] = mean(simunmax^2*(simun[,1] - simun[,2])^2)*K^2/4;
      varsimu2[1] = mean((simun[,1]*simun[,2]*4 + simun[,3]*(simun[,1] + simun[,2]))*(simusum+1)/12);
      for(ss in 2:(K-1)){
        simuss = apply(simun[,1:ss], 1, sum);
        varsimu1[ss] = mean(simunmax^2*(simuss - simun[,ss+1]*ss)^2)*K^2/4;
        varsimu2[ss] = mean((simuss*simun[, ss+1]*(ss+1)^2 + (simusum - simuss - simun[,ss+1])*(simuss + ss^2*simun[,ss+1]))*(simusum+1)/12);
      }
      varsimu = varsimu1 + varsimu2;
      w[jj, ii, kk] = sum(u^2/varsimu);
    }
    pw = pchisq(w, 2, lower.tail = F);
    powerkw[jj, ii] = mean(pkw[jj, ii, ] < 0.001);
    powerw[jj, ii] = mean(pw[jj, ii, ] < 0.001);
  }
  print(ii)
}


save(pw, pkw, powerkw, powerw, alphamatrix,N, betamatrix, theta, K, rep, file = "ANOVAexp2_seq.RData")


#################################################################################
#Exp 2c-d, equal size tests with the same nonzero functions and different theta#
##############################################################################
rep = 10000;
K = 3;

thetamatrix = c(0.1, 0.15, 0.22);
alpha = 1/2; beta = 2;
pkw = array(0, dim = c(2, 15, rep));
pw = array(0, dim = c(2, 15, rep));
powerkw = matrix(0, nrow = 2, ncol = 15);
powerw = matrix(0, nrow = 2, ncol = 15);
kw = pkw; w = pw;

for(ii in 1:15){
  N = ii*20;
  for(jj in 1:2){
    theta = thetamatrix + 0.05*(jj-1);
    for(kk in 1:rep){
      set.seed(20220*ii+10000*jj+kk)
      ###########Generate Data ###################
      X = rbeta(3*N, alpha, beta);
      Q1 = rbinom(N, 1, theta[1]);
      X1 = X[1:N]*Q1;
      Q2 = rbinom(N, 1, theta[2]);
      X2 = X[(N + 1):(2*N)]*Q2;
      Q3 = rbinom(N, 1, theta[3]);
      X3 = X[(2*N + 1):(3*N)]*Q3;
      
      ###some reads are missing ####
      ##when the reads are small, there is 1/2 probability missing ##
      ## small means the value is smaller than a ###
      a = 0.001;
      M = rbinom(length(X1), 1, 0.5);
      M[X1 > a] = 1;
      X1 = X1*M; 
      Q1 = 1*(X1 != 0);
      
      M = rbinom(length(X2), 1, 0.5);
      M[X2 > a] = 1;
      X2 = X2*M; 
      Q2 = 1*(X2 != 0);
      
      M = rbinom(length(X3), 1, 0.5);
      M[X3 > a] = 1;
      X3 = X3*M; 
      Q3 = 1*(X3 != 0);
      
      ###########Kruskal Wallis Test ############
      data = list(X1, X2, X3);
      pkw[jj, ii, kk] = kruskal.test(data)$p.value;
      
      ###########Our Test ###################
      #      pw[jj, ii, kk] = KW.zeros(data, g = NULL, alpha = 0.05, perm = TRUE)$p.value;
      n = c(sum(Q1), sum(Q2), sum(Q3)); nmax = max(n);
      X1trun = c(X1[Q1 == 1], rep(0, nmax - n[1]));
      X2trun = c(X2[Q2 == 1], rep(0, nmax - n[2]));
      X3trun = c(X3[Q3 == 1], rep(0, nmax - n[3]));
      rankdata = K*nmax + 1 - rank(c(X1trun, X2trun, X3trun));
      r = c(sum(rankdata[1:nmax]), sum(rankdata[nmax+1:nmax]), sum(rankdata[2*nmax+1:nmax]) );
      s = r - nmax*(K*nmax + 1)/2;
      u = c(s[1] - s[2], sum(s[1:2]) - 2*s[3]);
      thetam = mean(n)/N;
      
      simun = matrix(rbinom(5000*K, N, thetam), nrow = 5000, ncol = K);
      simunmax = apply(simun, 1, max); simusum = apply(simun, 1, sum);
      varsimu1 = numeric(K - 1); varsimu2 = varsimu1;
      varsimu1[1] = mean(simunmax^2*(simun[,1] - simun[,2])^2)*K^2/4;
      varsimu2[1] = mean((simun[,1]*simun[,2]*4 + simun[,3]*(simun[,1] + simun[,2]))*(simusum+1)/12);
      for(ss in 2:(K-1)){
        simuss = apply(simun[,1:ss], 1, sum);
        varsimu1[ss] = mean(simunmax^2*(simuss - simun[,ss+1]*ss)^2)*K^2/4;
        varsimu2[ss] = mean((simuss*simun[, ss+1]*(ss+1)^2 + (simusum - simuss - simun[,ss+1])*(simuss + ss^2*simun[,ss+1]))*(simusum+1)/12);
      }
      varsimu = varsimu1 + varsimu2;
      w[jj, ii, kk] = sum(u^2/varsimu);
    }
    pw = pchisq(w, 2, lower.tail = F);
    powerkw[jj, ii] = mean(pkw[jj, ii, ] < 0.001);
    powerw[jj, ii] = mean(pw[jj, ii, ] < 0.001);
  }
  print(ii)
}

save(w, pkw, pw, powerkw, powerw, alpha, beta,N, thetamatrix, rep, file = "ANOVAexp3_seq.RData")
rm(list = ls())





#############################################################################
#Exp 2e-f, unequal size tests with the different nonzero functions and different theta
#############################################################################

rep = 10000;
K = 3;

theta = c(0.4, 0.5, 0.6); Nprop = c(0.8, 1, 1.5);
alphamatrix = matrix(c(1.5, 2, 2.5, 1.7, 2, 2.3, 1.5, 2, 2.5, 1.7, 2, 2.3), nrow = 4, ncol = K, byrow = T);
betamatrix = matrix(c(2, 2, 2, 2, 2, 2, 4, 4, 4, 4, 4, 4), nrow = 4, ncol = K, byrow = T);
pkw = array(0, dim = c(2, 20, rep));
pw = array(0, dim = c(2, 20, rep));
powerkw = matrix(0, nrow = 2, ncol = 20);
powerw = matrix(0, nrow = 2, ncol = 20);
w = pw;

for(ii in 1:20){
  Nbase = ii*100;
  for(jj in 1:2){
    alpha = alphamatrix[jj,];
    beta = betamatrix[jj, ];
    for(kk in 1:rep){
      set.seed(20220*ii+10000*jj+kk)
      N = round(Nbase*Nprop);
      ###########Generate Data ###################
      Q = rbinom(sum(N), 1, mean(theta));
      X1 = rbeta(N[1], alpha[1], beta[1])
      X1 = Q[1:N[1]]*X1;
      X2 = rbeta(N[2], alpha[2], beta[2])
      X2 = Q[N[1] + 1:N[2]]*X2;
      X3 = rbeta(N[3], alpha[3], beta[3])
      X3 = Q[N[1] + N[2] +1:N[3]]*X3;
      
      ###some reads are missing ####
      ##when the reads are small, there is 1/2 probability missing ##
      ## small means the value is smaller than a ###
      a = 0.1;
      M = rbinom(length(X1), 1, 0.5);
      M[X1 > a] = 1;
      X1 = X1*M; 
      Q1 = 1*(X1 != 0);
      
      M = rbinom(length(X2), 1, 0.5);
      M[X2 > a] = 1;
      X2 = X2*M; 
      Q2 = 1*(X2 != 0);
      
      M = rbinom(length(X3), 1, 0.5);
      M[X3 > a] = 1;
      X3 = X3*M; 
      Q3 = 1*(X3 != 0);
      
      Q = c(Q1, Q2, Q3);
      ###########Kruskal Wallis Test ############
      data = list(X1, X2, X3);
      pkw[jj, ii, kk] = kruskal.test(data)$p.value;
      
      ###########Our Test ###################
      n = c(sum(Q[1:N[1]]), sum(Q[N[1] + 1:N[2]]), sum(Q[N[1] + N[2] +1:N[3]]))
      prop = n/N; pmax = max(prop);
      Ntrun = round(pmax*N);
      X1trun = c(X1[Q[1:N[1]] == 1], rep(0, Ntrun[1] - n[1]));
      X2trun = c(X2[Q[N[1] + 1:N[2]] == 1], rep(0, Ntrun[2]- n[2]));
      X3trun = c(X3[Q[N[1] + N[2] +1:N[3]] == 1], rep(0, Ntrun[3]- n[3]));
      rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun, X3trun));
      r = c(sum(rankdata[1:Ntrun[1]]), sum(rankdata[Ntrun[1]+1:Ntrun[2]]), sum(rankdata[sum(Ntrun[1:2])+1:Ntrun[3]]) ); 
      s = r - Ntrun*(sum(Ntrun) + 1)/2;
      u = c(N[2]*s[1] - N[1]*s[2], N[3]*sum(s[1:2]) - s[3]*sum(N[1:2]))/Nbase^2;
      
      thetam = mean(prop);
      simun = matrix(0, nrow = 5000, ncol = K); simup = simun;
      for(ss in 1:K){
        simun[,ss] = rbinom(5000, N[ss], thetam);		
        simup[,ss] = simun[,ss]/N[ss];
      }
      simupmax = apply(simup, 1, max);
      varsimu = numeric(K - 1);
      varsimu[1] = Nprop[2]^2*mean(simupmax^2*(simup[,1] - simup[,2])^2)*N[1]^2;
      for(ss in 2:(K-1)){
        varsimu[ss] = Nprop[ss+1]^2*mean(simupmax^2*(apply(simun[,1:ss], 1, sum) - simup[,ss+1]*sum(N[1:ss]))^2);
      }
      varsimu = varsimu*(sum(Nprop))^2/4;
      
      varu2 = c(Nprop[1]*Nprop[2], Nprop[3]*sum(Nprop))*(Nprop[1] + Nprop[2])*sum(Nprop)*thetam^2*(sum(N)*thetam + 1)/12;
      varu = varsimu + varu2;
      w[jj, ii, kk] = sum(u^2/varu);
    }
    pw = pchisq(w, 2, lower.tail = F);
    powerkw[jj, ii] = mean(pkw[jj, ii, ] < 0.001);
    powerw[jj, ii] = mean(pw[jj, ii, ] < 0.001);
  }
  print(ii)
}
save(pw, pkw, powerkw, powerw, alphamatrix,N, Nprop, betamatrix, theta, K, rep, file = "ANOVAexp4_seq.RData")



### Draw Plots ####

pdf("power.pdf",width=6,height=9)

par(mfrow = c(3,2))

for(jj in 1:2){
  load("ANOVAexp2_seq.Rdata")
  N = (1:20)*100;
  plot(N, powerkw[jj,], type = "l", lwd = 2, lty = 2, ylim = c(0, 1), ylab = "Power", col = 2)
  lines(N, powerw[jj,], lwd = 2, col = 2)
  load("ANOVAexp2.Rdata")
  N = (1:20)*100;
  lines(N, powerkw[jj,], lwd = 2, lty = 2)
  lines(N, powerw[jj,], lwd = 2)
  if(jj == 1){text(200, 0.9, "(a)", pos = 3, cex = 1.2)}
  if(jj == 2){text(200, 0.9, "(b)", pos = 3, cex = 1.2)}
}

for(jj in 1:2){
  load("ANOVAexp3_seq.Rdata")
  N = (1:15)*20;
  plot(N, powerkw[jj,1:15], type = "l", lwd = 2, lty = 2, ylim = c(0, 1), ylab = "Power", col = 2)
  lines(N, powerw[jj,1:15], lwd = 2, col = 2)
  load("ANOVAexp3.Rdata")
  N = (1:15)*20;
  lines(N, powerkw[jj,1:15], lwd = 2, lty = 2)
  lines(N, powerw[jj,1:15], lwd = 2)
  if(jj == 1){text(30, 0.9, "(c)", pos = 3, cex = 1.2)}
  if(jj == 2){text(30, 0.9, "(d)", pos = 3, cex = 1.2)}
}

for(jj in 1:2){
  load("ANOVAexp4_seq.Rdata")
  N = (1:20)*100;
  plot(N, powerkw[jj,], type = "l", lwd = 2, lty = 2, ylim = c(0, 1), ylab = "Power", col = 2)
  lines(N, powerw[jj,], lwd = 2, col = 2)
  load("ANOVAexp4.Rdata")
  N = (1:20)*100;
  lines(N, powerkw[jj,], lwd = 2, lty = 2)
  lines(N, powerw[jj,], lwd = 2)
  if(jj == 1){text(200, 0.9, "(e)", pos = 3, cex = 1.2)}
  if(jj == 2){text(200, 0.9, "(f)", pos = 3, cex = 1.2)}
}
dev.off()
