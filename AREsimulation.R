### This file is the simulation file to verify the ARE curve for truncated 
## Wilcoxon statistic and standard statistic ##############################

#############################################################################
#setting (a). Effect of Delta_{f,g}##########################################
#############################################################################

rep = 10000;

theta1 = 0.3; theta2 = 0.8;
#for each delta_{f,g}, we want to find proper Beta dist parameters
#fix alpha1 and beta1, change beta2
alpha1 = 3; beta1 = 3; 
alpha.table = seq(1, 10, by = 0.01)
beta.table = seq(1, 10, by = 0.01)
delta.table = matrix(0, nrow = length(alpha.table), ncol = length(beta.table)); 

xx = seq(0, 1, by = 0.001); 
for(ii in 1:length(alpha.table)){
  for(jj in 1:length(beta.table)){
    delta.table[ii, jj] = sum(dbeta(xx, alpha.table[ii], beta.table[jj])*
                                pbeta(xx, alpha1, beta1, lower.tail = TRUE)*0.001) - 1/2;
  }
}


N = c(40, 50);
deltaarray = seq(-0.47, 0.48, by = 0.01)
kw = matrix(0, nrow = length(deltaarray), ncol = rep);
pkw = kw; skw = kw; varkw = kw; 
w = matrix(0, nrow = length(deltaarray), ncol = rep);
sw = w; varw = w; 
rekw = rep(0, length(deltaarray));
rew = rep(0, length(deltaarray));

for(ii in 1:length(deltaarray)){
  delta = deltaarray[ii]; 
  result = which(abs(delta - delta.table) == min(abs(delta - delta.table)), arr.ind = TRUE);
  alpha2 = alpha.table[result[1]];
  beta2 = beta.table[result[2]];
  for(kk in 1:rep){
    ###########Generate Data ###################
    X1 = rbeta(N[1], alpha1, beta1)
    Q1 = rbinom(N[1], 1, theta1); 
    X1 = Q1*X1;
    X2 = rbeta(N[2], alpha2, beta2)
    Q2 = rbinom(N[2], 1, theta2); 
    X2 = Q2*X2;
    
    ###########Wilcoxon Rank Sum Test ############
    pkw[ii, kk] = wilcox.test(X1, X2, correct = FALSE)$p.value
    kw[ii, kk] = qnorm(1 - pkw[ii, kk]/2, lower.tail = TRUE)
    
    ###########Our Test ###################
    n = c(sum(Q1), sum(Q2))
    prop = n/N; pmax = max(prop); pmean = mean(prop);
    Ntrun = round(pmax*N);
    X1trun = c(X1[Q1 == 1], rep(0, Ntrun[1] - n[1]));
    X2trun = c(X2[Q2 == 1], rep(0, Ntrun[2]- n[2]));
    rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun));
    r = sum(rankdata[1:Ntrun[1]]);
    s = r - Ntrun[1]*(sum(Ntrun) + 1)/2 - pmean*(1 - pmean)*(N[2]-N[1])/4;
    sw[ii, kk] = s;
    
    #### Variance part ##############
    ### Truncated Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*pmean*v1*sum(N)*pmean - v1^2*(N[1] - N[2])^2/4;
    var1 = var1/4;
    var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
    vars = var1 + var2;
    varw[ii, kk] = vars;
    
    ### Original Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*v1*sum(N); 
    var1 = var1/4;
    varS = var1 + var2;
    varkw[ii, kk] = varS;
    
    ### Test Statistic ###
    w[ii, kk] = s/sqrt(vars);
  }
  print(ii)
}

for(ii in 1:length(deltaarray)){
  rekw[ii] = mean(kw[ii, ]^2);
  rew[ii] = mean(w[ii, ]^2);
}

plot(deltaarray, rew/rekw, type = "l", lwd = 2, xlab = expression(Delta*theta), ylab = expression("ARE(T"[tW]*", T"[W]*")"))
abline(h = 1, lty = 2)
text(0.05, 2.45, "(a)", pos = 3, cex = 1.2)

save(w, kw, rew, rekw, deltaarray, alpha1, beta1, alpha2, beta2, rep, N, file = "AREsimua.RData")
rm(list = ls())



#############################################################################
#setting (b). Effect of theta##########################################
#############################################################################

rep = 10000;

theta.all = seq(0.1, 0.9, by = 0.01);
alpha = 2; beta1 = 2.75; beta2 = 2;
#show that delta_{f,g} = 0.1
xx = seq(0, 1, by = 0.001); 
sum(dbeta(xx, alpha, beta2)*pbeta(xx, alpha, beta1, lower.tail = TRUE)*0.001)

N = c(40, 50);
kw = matrix(0, nrow = length(theta.all), ncol = rep);
pkw = kw; skw = kw; varkw = kw; 
w = matrix(0, nrow = length(theta.all), ncol = rep);
sw = w; varw = w; 
rekw = rep(0, length(theta.all));
rew = rep(0, length(theta.all));

for(ii in 1:length(theta.all)){
  theta = theta.all[ii]; 
  theta1 = theta - 0.1; theta2 = theta + 0.1;
  for(kk in 1:rep){
    ###########Generate Data ###################
    X1 = rbeta(N[1], alpha, beta2)
    Q1 = rbinom(N[1], 1, theta1); 
    X1 = Q1*X1;
    X2 = rbeta(N[2], alpha, beta1)
    Q2 = rbinom(N[2], 1, theta2); 
    X2 = Q2*X2;
    
    ###########Wilcoxon Rank Sum Test ############
    pkw[ii, kk] = wilcox.test(X1, X2, correct = FALSE)$p.value
    kw[ii, kk] = qnorm(1 - pkw[ii, kk]/2, lower.tail = TRUE)
    
    ###########Our Test ###################
    n = c(sum(Q1), sum(Q2))
    prop = n/N; pmax = max(prop); pmean = mean(prop);
    Ntrun = round(pmax*N);
    X1trun = c(X1[Q1 == 1], rep(0, Ntrun[1] - n[1]));
    X2trun = c(X2[Q2 == 1], rep(0, Ntrun[2]- n[2]));
    rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun));
    r = sum(rankdata[1:Ntrun[1]]);
    s = r - Ntrun[1]*(sum(Ntrun) + 1)/2 - pmean*(1 - pmean)*(N[2]-N[1])/4;
    sw[ii, kk] = s;
    
    #### Variance part ##############
    ### Truncated Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*pmean*v1*sum(N)*pmean - v1^2*(N[1] - N[2])^2/4;
    var1 = var1/4;
    var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
    vars = var1 + var2;
    varw[ii, kk] = vars;
    
    ### Original Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*v1*sum(N); 
    var1 = var1/4;
    varS = var1 + var2;
    varkw[ii, kk] = varS;
    
    ### Test Statistic ###
    w[ii, kk] = s/sqrt(vars);
#    kw[ii, kk] = S/sqrt(varS); 
  }
  print(ii)
}

for(ii in 1:length(theta.all)){
  rekw[ii] = mean(kw[ii, ]^2);
  rew[ii] = mean(w[ii, ]^2);
}

 save(w, kw, rew, rekw, theta.all, alpha, beta1, beta2, rep, N, file = "AREsimub.RData")
 rm(list = ls())


#############################################################################
#setting (c). Effect of Delta theta##########################################
#############################################################################

rep = 10000;

theta.delta = seq(0.02, 0.9, by = 0.02);
alpha = 2; beta1 = 2.75; beta2 = 2;
#show that delta_{f,g} = 0.1
xx = seq(0, 1, by = 0.001); 
sum(dbeta(xx, alpha, beta2)*pbeta(xx, alpha, beta1, lower.tail = FALSE)*0.001)

N = c(40, 50);
kw = matrix(0, nrow = length(theta.delta), ncol = rep);
skw = kw; varkw = kw; pkw = kw;
w = matrix(0, nrow = length(theta.delta), ncol = rep);
sw = w; varw = w; 
rekw = rep(0, length(theta.delta));
rew = rep(0, length(theta.delta));


for(ii in 1:length(theta.delta)){
  dtheta = theta.delta[ii]; 
  theta1 = 0.5 - dtheta/2; theta2 = 0.5 + dtheta/2;
  for(kk in 1:rep){
    ###########Generate Data ###################
    X1 = rbeta(N[1], alpha, beta1)
    Q1 = rbinom(N[1], 1, theta1); 
    X1 = Q1*X1;
    X2 = rbeta(N[2], alpha, beta2)
    Q2 = rbinom(N[2], 1, theta2); 
    X2 = Q2*X2;
    
    ###########Wilcoxon Rank Sum Test ############
    pkw[ii, kk] = wilcox.test(X1, X2, correct = FALSE)$p.value
    kw[ii, kk] = qnorm(1 - pkw[ii, kk]/2, lower.tail = TRUE)
    if(pkw[ii, kk] < 1e-15) {
      rank1 = sum(N) + 1 - rank(c(X1, X2));
      R = sum(rank1[1:N[1]]);
      S = R - N[1]*(sum(N) + 1)/2;
      skw[ii, kk] = S;
    }
    
    
    ###########Our Test ###################
    n = c(sum(Q1), sum(Q2))
    prop = n/N; pmax = max(prop); pmean = mean(prop);
    Ntrun = round(pmax*N);
    X1trun = c(X1[Q1 == 1], rep(0, Ntrun[1] - n[1]));
    X2trun = c(X2[Q2 == 1], rep(0, Ntrun[2]- n[2]));
    rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun));
    r = sum(rankdata[1:Ntrun[1]]);
    s = r - Ntrun[1]*(sum(Ntrun) + 1)/2 - pmean*(1 - pmean)*(N[2]-N[1])/4;
    sw[ii, kk] = s;
    
    #### Variance part ##############
    ### Truncated Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*pmean*v1*sum(N)*pmean - v1^2*(N[1] - N[2])^2/4;
    var1 = var1/4;
    var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
    vars = var1 + var2;
    varw[ii, kk] = vars;
    
    ### Original Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*v1*sum(N); 
    var1 = var1/4;
    varS = var1 + var2;
    varkw[ii, kk] = varS;
    
    ### Test Statistic ###
    w[ii, kk] = s/sqrt(vars);
    if(pkw[ii, kk] < 1e-15) {kw[ii, kk] = S/sqrt(varS);} 
  }
  print(ii)
}
for(ii in 1:length(theta.delta)){
  rekw[ii] = mean(kw[ii, ]^2);
  rew[ii] = mean(w[ii, ]^2);
}

save(w, kw, rew, rekw, theta.delta, alpha, beta1, beta2, rep, N, file = "AREsimuc.RData")
rm(list = ls())

#############################################################################
#setting (d). Effect of N1##########################################
#############################################################################

rep = 10000;

Narray = seq(20, 120, by = 5);
theta1 = 0.3; theta2 = 0.8;
alpha = 2; beta1 = 2.75; beta2 = 2;
#show that delta_{f,g} = 0.1
xx = seq(0, 1, by = 0.001); 
sum(dbeta(xx, alpha, beta2)*pbeta(xx, alpha, beta1, lower.tail = TRUE)*0.001)

kw = matrix(0, nrow = length(Narray), ncol = rep);
skw = kw; varkw = kw; pkw = kw;
w = matrix(0, nrow = length(Narray), ncol = rep);
sw = w; varw = w; 

for(ii in 1:length(Narray)){
  N = c(Narray[ii], 50);
  
  for(kk in 1:rep){
    ###########Generate Data ###################
    X1 = rbeta(N[1], alpha, beta1)
    Q1 = rbinom(N[1], 1, theta1); 
    X1 = Q1*X1;
    X2 = rbeta(N[2], alpha, beta2)
    Q2 = rbinom(N[2], 1, theta2); 
    X2 = Q2*X2;
    
    ###########Wilcoxon Rank Sum Test ############
    pkw[ii, kk] = wilcox.test(X1, X2, correct = FALSE)$p.value
    kw[ii, kk] = qnorm(1 - pkw[ii, kk]/2, lower.tail = TRUE)
    if(pkw[ii, kk] < 1e-15) {
      rank1 = sum(N) + 1 - rank(c(X1, X2));
      R = sum(rank1[1:N[1]]);
      S = R - N[1]*(sum(N) + 1)/2;
      skw[ii, kk] = S;
    }
    
    
    ###########Our Test ###################
    n = c(sum(Q1), sum(Q2))
    prop = n/N; pmax = max(prop); pmean = mean(prop);
    Ntrun = round(pmax*N);
    X1trun = c(X1[Q1 == 1], rep(0, Ntrun[1] - n[1]));
    X2trun = c(X2[Q2 == 1], rep(0, Ntrun[2]- n[2]));
    rankdata = sum(Ntrun) + 1 - rank(c(X1trun, X2trun));
    r = sum(rankdata[1:Ntrun[1]]);
    s = r - Ntrun[1]*(sum(Ntrun) + 1)/2 - pmean*(1 - pmean)*(N[2]-N[1])/4;
    sw[ii, kk] = s;
    
    #### Variance part ##############
    ### Truncated Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*pmean*v1*sum(N)*pmean - v1^2*(N[1] - N[2])^2/4;
    var1 = var1/4;
    var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
    vars = var1 + var2;
    varw[ii, kk] = vars;
    
    ### Original Test ####
    v1 = pmean*(1 - pmean);
    var1 = N[1]*N[2]*v1*sum(N); 
    var1 = var1/4;
    varS = var1 + var2;
    varkw[ii, kk] = varS;
    
    ### Test Statistic ###
    w[ii, kk] = s/sqrt(vars);
    if(pkw[ii, kk] < 1e-15) {kw[ii, kk] = S/sqrt(varS);} 
  }
  print(ii)
}

rekw = apply(kw^2, 1, mean);
rew = apply(w^2, 1, mean);

plot(Narray, rew/rekw, type = "l", lwd = 2, xlab = expression(Delta*theta), ylab = expression("ARE(T"[tW]*", T"[W]*")"))

save(w, kw, rew, rekw, theta1, theta2, alpha, beta1, beta2, rep, Narray, file = "AREsimud.RData")
rm(list = ls())



### Draw Overlay Plot ####
pdf(file = "ARE_simulation.pdf", width = 8, height = 8)
par(mfrow = c(2, 2))

N1 = 40; N2 = 50;
theta1 = 0.3; theta2 = 0.8;
n1 = N1*theta1; n2 = N2*theta2;
thetam = max(theta1, theta2);
theta = (theta1 + theta2)/2;
dtheta = theta2 - theta1;
d12 = seq(-0.5, 0.5, 0.01);

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

par(mar=c(4,4,2,2)) 
load("AREsimua.RData")
plot(d12, ratio, type = "l", lwd = 2, ylim = c(min(ratio, rew/rekw), max(ratio, rew/rekw)), 
     xlab = expression(Delta["f,g"]), 
     ylab = expression("ARE(T"[tW]*", T"[W]*")"))
text(x = -0.48, y = 1.7, "(a)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
# text(-0.48, 1.05, "y=1", cex = 1)
lines(deltaarray, rew/rekw, lwd = 2, col = 2)


N1 = 40; N2 = 50;
theta = seq(0.1, 0.9, by = 0.01);
theta1 = theta - 0.1; theta2 = theta + 0.1;
n1 = N1*theta1; n2 = N2*theta2;
thetam = apply(cbind(theta1, theta2), 1, max);
theta = (theta1 + theta2)/2;
dtheta = theta2 - theta1;
d12 = 0.1;

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

load("AREsimub.RData")
par(mar=c(4,4,2,4)) 
plot(theta, ratio, type = "l", lwd = 2, ylim = c(min(ratio, rew/rekw), max(ratio, rew/rekw)), 
     xlab = expression(theta), ylab = expression("ARE(T"[tW]*", T"[W]*")"))
text(0.85, 2.95, "(b)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
# text(0.15, 1.05, "y=1", cex = 1)
lines(theta.all, rew/rekw, lwd = 2, col = 2)


N1 = 40; N2 = 50;
dtheta = seq(0, 0.5, by = 0.01);
theta1 = 0.5 - dtheta; theta2 = 0.5 + dtheta;
n1 = N1*theta1; n2 = N2*theta2;
thetam = apply(cbind(theta1, theta2), 1, max);
theta = (theta1 + theta2)/2;
d12 = 0.1;

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);


load("AREsimuc.RData")
par(mar=c(4,4,2,2)) 
plot(dtheta*2, ratio, type = "l", lwd = 2, ylim = c(min(ratio, rew/rekw), max(ratio, rew/rekw)), 
     xlab = expression(Delta*theta), ylab = expression("ARE(T"[tW]*", T"[W]*")"))
text(0.03, 2.65, "(c)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
lines(theta.delta, rew/rekw, lwd = 2, col = 2)

N1 = seq(20, 120, by = 5); N2 = 50;

theta1 = 0.3; theta2 = 0.8;
n1 = N1*theta1; n2 = N2*theta2;
thetam = max(theta1, theta2);
theta = (theta1 + theta2)/2;
dtheta = theta1 - theta2;
d12 = 0.1;

vS = N1*N2*(N1+N2)*theta*(theta^2+3*(1 - theta))/12;
vs = N1*N2*(N1+N2)*theta^3*(1+3*(1 - theta))/12;
meanS = n1*n2*d12 + N1*N2*dtheta/2;
means = n1*n2*d12 + N1*N2*dtheta/2*thetam - theta*(1 - theta)*(N2 - N1)/4;
ratio = (means/meanS)^2*(vS/vs);

ratio = (1/theta^2)*(1 - theta + theta^2/3)/(1 - theta + 1/3)*((theta1*theta2*d12 + dtheta*thetam/2)/(theta1*theta2*d12 + dtheta/2))^2
load("AREsimud.RData")

par(mar=c(4,4,2,4))
plot(N1, rep(ratio, length(N1)), type = "l", lwd = 2, ylim = c(min(ratio, rew/rekw), max(ratio, rew/rekw)),
     xlab = expression('N'[1]), ylab = expression("ARE(T"[tW]*", T"[W]*")"))
text(116, 1.75, "(d)", pos = 3, cex = 1.2)
abline(h = 1, lty = 2)
lines(Narray, rew/rekw, lwd = 2, col = 2)
dev.off()