#########################################################
# netW.R
#
# The version of the Wagner model used in 
# Rünneburger & Le Rouzic 2016 BMC Evol. Biol. 
#
# 
# Copyright Arnaud Le Rouzic / CNRS 2015-2017
# <lerouzic@egce.cnrs-gif.fr>
#
# Released under the WTFPL version 2.0
# * No warranty *
###########################################################

sigma.M2 <- function(x, a) {
	1. / (1. + exp((-x/(a*(a-1)))+log(1/a-1)))
}

sigma.M2p <- function(x, lambda, mu) {
	1. / (1. + lambda * exp(-mu*x))
}

suppressMessages(library(compiler))
sigma.M2c <- cmpfun(sigma.M2p)

#~ library(Rcpp)

#~ cppFunction('
#~ NumericVector sigma_M2p(NumericVector x, double aam1, double l1am1) {
#~     NumericVector ans(x.size());
#~     // double aam1 = a*(1.-a);
#~     // double l1am1 = log(1./a-1.);
#~     for (int i = 0; i < x.size(); i++) {
#~         ans[i] = 1. / (1. + exp((-x[i]/aam1)+l1am1));
#~     }
#~     return ans;
#~ }
#~ ')

#~ cppFunction('
#~ NumericVector sigma_M2(NumericVector x, double a) {
#~     NumericVector ans(x.size());
#~     double aam1 = a*(1.-a);
#~     double l1am1 = log(1./a-1.);
#~     for (int i = 0; i < x.size(); i++) {
#~         ans[i] = 1. / (1. + exp((-x[i]/aam1)+l1am1));
#~     }
#~     return ans;
#~ }
#~ ')

internal_loop_R <- function(W, S0, a, steps, measure) {
	lambda <- (1-a)/a
	mu <- 1/(a*(1-a))
	sto <- matrix(NA, nrow=length(S0), ncol=steps+1)
	sto[,1] <- S0
	for (i in 1:steps) {
		S0 <- sigma.M2c((W %*% S0), lambda=lambda, mu=mu) 			
		sto[,i+1] <- S0
	}
	list(full=sto, sumx=rowSums(sto[,(steps-measure+1):steps]), sumx2=rowSums(sto[,(steps-measure+1):steps]^2))
}

library(Rcpp)
library(inline, quietly=TRUE)

cppFunction('
	List internal_loop_cpp(const NumericMatrix &W, const NumericVector &S0, double a, unsigned int steps, unsigned int measure) {
		double lambda = (1-a)/a;
		double mu     = 1/(a*(1-a));
		NumericMatrix sto (S0.size(), steps+1);
		NumericVector sumx (S0.size());
		NumericVector sumx2 (S0.size()); 
		for (unsigned int i = 0; i < S0.size(); i++)
			sto(i,0) = S0(i);
		for (unsigned int t = 1; t <= steps; t++) {
			for (unsigned int i = 0; i < S0.size(); i++) {
				double tmp = 0.;
				for (unsigned int j = 0; j < S0.size(); j++) {
					tmp += sto(j,t-1) * W(i,j);
				}
				tmp =  1. / (1. + lambda * exp(-mu*tmp));
				sto(i,t) = tmp;
				if (t > steps-measure) {
					sumx(i) += tmp;
					sumx2(i) += tmp*tmp;
				}
			}
		}
		for (unsigned int i = 0; i < S0.size(); i++) {
			sumx(i) /= static_cast<double>(measure); // sumx(i) now contains the mean
			sumx2(i) /= static_cast<double>(measure);
			sumx2(i) += -sumx(i)*sumx(i); // sumx2(i) now contains the variance
		}
		return List::create(Named("full")=sto, Named("mean")=sumx, Named("var")=sumx2);
	}')


model.M2 <- function(W, a=0.5, S0=rep(a, nrow(W)), steps=20, measure=4, full=FALSE, loopFUN=internal_loop_cpp) {
    ans <- loopFUN(W, S0, a, steps, measure)
	if (!full) ans$full <- NULL
	return(ans)
}

model.M2 <- cmpfun(model.M2)
