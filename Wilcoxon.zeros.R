# This function is to calculate the Wilcoxon Rak Sum statistic when the data has 
# a clump of zeros. 
#

# Inputs:
# x: a numeric vector of data values from one sample.
# y: a numeric vector of data values from another sample.
# alpha: a constant giving the level of the testing problem. Default value as 0.05
# perm: TRUE/FALSE value. If TURE, then the p-value will be calculated through permutations. Default as FALSE.
#
# Outputs:
# H: TRUE/FALSE value, giving the hypothesis testing result
# p.value: the p-value 
# statistics: the modified Wilcoxon Rank Sum statistic

Wilcoxon.zeros <- function(x, y, alpha = 0.05, perm = FALSE, numrep = 10000) {
	
	N = c(length(x), length(y));

	###########Our Test ###################
	n = c(sum(x != 0), sum(y != 0));
	prop = n/N; pmax = max(prop); pmean = mean(prop);
	Ntrun = round(pmax*N);
	Xtrun = c(x[x != 0], rep(0, Ntrun[1] - n[1]));
	Ytrun = c(y[y != 0], rep(0, Ntrun[2]- n[2]));
	rankdata = sum(Ntrun) + 1 - rank(c(Xtrun, Ytrun));
	r = sum(rankdata[1:Ntrun[1]]);
	s = r - Ntrun[1]*(sum(Ntrun) + 1)/2 - pmean*(1 - pmean)*(N[2]-N[1])/4;

	v1 = pmean*(1 - pmean);
	# var1 = N[1]*N[2]*pmean*v1*(sum(N)*pmean + 3*(1 - pmean)/2) + sum(N^2)*v1^2*5/4 + 2*sqrt(2/pi)*pmean*(sum(N)*v1)^(1.5)*sqrt(N[1]*N[2]);
	var1 = N[1]*N[2]*pmean*v1*sum(N)*pmean - v1^2*(N[1] - N[2])^2/4;
	var1 = var1/4;
	var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
	vars = var1 + var2;
	w = s/sqrt(vars);

	if(perm == FALSE) {p = 2*(1 - pnorm(abs(w)));}
	else{
		permu.w = rep(0, numrep); 
		Z = c(x, y);
		for(i in 1:numrep){
			ind = sample(sum(N), N[1]);
			x = Z[ind]; y = Z[-ind];
			n = c(sum(x != 0), sum(y != 0));
			prop = n/N; pmax = max(prop); pmean = mean(prop);
			Ntrun = round(pmax*N);
			Xtrun = c(x[x != 0], rep(0, Ntrun[1] - n[1]));
			Ytrun = c(y[y != 0], rep(0, Ntrun[2]- n[2]));
			rankdata = sum(Ntrun) + 1 - rank(c(Xtrun, Ytrun));
			r = sum(rankdata[1:Ntrun[1]]);
			s = r - Ntrun[1]*(sum(Ntrun) + 1)/2 - pmean*(1 - pmean)*(N[2]-N[1])/4;

			v1 = pmean*(1 - pmean);
#			var1 = N[1]*N[2]*pmean*v1*(sum(N)*pmean + 3*(1 - pmean)/2) + sum(N^2)*v1^2*5/4 + 2*sqrt(2/pi)*pmean*(sum(N)*v1)^(1.5)*sqrt(N[1]*N[2]);
			var1 = N[1]*N[2]*pmean*v1*sum(N)*pmean - v1^2*(N[1] - N[2])^2/4;
			var1 = var1/4;
			var2 = N[1]*N[2]*pmean^2*(sum(N)*pmean + 1)/12;
			vars = var1 + var2;
			permu.w[i] = s/sqrt(vars);
		}
		p = sum(abs(w) < abs(permu.w))/numrep;	
	}
	H = (p < alpha);

	result <- list(H = H, p.value = p, statistics=w^2)
	return(result)
}
