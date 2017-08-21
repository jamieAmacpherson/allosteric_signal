#! /usr/bin/R

#===============================================================================
# compare MI matrices
#===============================================================================

library("entropy");

#_______________________________________________________________________________
## 5x5 toy system
## toy matrix 'a'
a.t.m = c(0.0, 0.9, 0.0, 0.5, 0.0,
		0.9, 0.0, 0.0, 0.0, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0,
		0.5, 0.0, 0.0, 0.0, 0.0,
		0.0, 0.5, 0.0, 0.0, 0.0);
dim(a.t.m) = c(5, 5);
image(a.t.m);

## toy matrix 'b'
## like 'a.m', but with two additional correlations in position 3
b.t.m = c(0.0, 0.9, 0.0, 0.5, 0.0,
		0.9, 0.0, 0.7, 0.0, 0.5,
		0.0, 0.7, 0.0, 0.7, 0.0,
		0.5, 0.0, 0.7, 0.0, 0.0,
		0.0, 0.5, 0.0, 0.0, 0.0);
dim(b.t.m) = c(5, 5);
image(b.t.m);

#_______________________________________________________________________________
## difference matrix
c.t.m = b.t.m - a.t.m;
image(c.t.m);

## Conclusion: For sa system comparison, spectral decomposition should be done
##   on the difference matrix, because equal signals in the matrices
##   (and the systems they represent cancel out and are not of interest.

## Of course the spectral decomposition of the single matrices is of interest
##   for the detection of correlated state transitions in systems 'a' and 'b'.

#_______________________________________________________________________________
#_______________________________________________________________________________
## 50x50 system with positive Gaussian noise
## insert MI values into given MI matrix
## 'mimat' is the MI matrix, 'posvec' the vector of positions to correlate and
##   'mival' the MI value to insert
insert_cor = function(mimat, posvec, mival) {
	cmb = combn(posvec, 2);
	for (k in 1:dim(cmb)[2]) {
		i = cmb[1, k];
		j = cmb[2, k];
		mimat[i, j] = mimat[i, j] + mival;
		mimat[j, i] = mimat[j, i] + mival;
	};
	return(mimat);
}

## Gaussian noise matrix 'a'
## create positive Gaussian noise
a.gn.m = abs(rnorm(2500));
## scale to [0,1]
a.gn.m = (a.gn.m - min(a.gn.m)) / (max(a.gn.m) - min(a.gn.m));
## set noise level
noise.level = 0.2;
a.gn.m = a.gn.m * noise.level;
dim(a.gn.m) = c(50, 50);
hist(a.gn.m);
## add correlations
## 0.1
c1 = c(2, 4, 6, 8, 10);
a.gn.m = insert_cor(a.gn.m, c1, 0.1);
## 0.3
c3 = c(1, 11, 21, 31, 41);
a.gn.m = insert_cor(a.gn.m, c3, 0.2);
## 0.5
c5 = c(12, 13, 14, 15, 16);
a.gn.m = insert_cor(a.gn.m, c5, 0.3);


## Gaussian noise matrix 'b'
## like 'a.gn.m', but with 3 fewer correlations in 'c1', 'c3' and 'c5'
## create positive Gaussian noise
b.gn.m = abs(rnorm(2500));
## scale to [0,1]
b.gn.m = (b.gn.m - min(b.gn.m)) / (max(b.gn.m) - min(b.gn.m));
## set noise level
noise.level = 0.2;
b.gn.m = b.gn.m * noise.level;
dim(b.gn.m) = c(50, 50);
hist(b.gn.m);
## add correlations
## 0.1
c1 = c(2, 10);
b.gn.m = insert_cor(b.gn.m, c1, 0.1);
## 0.3
c3 = c(1, 41);
b.gn.m = insert_cor(b.gn.m, c3, 0.2);
## 0.5
c5 = c(12, 16);
b.gn.m = insert_cor(b.gn.m, c5, 0.3);

#_______________________________________________________________________________
## difference matrix
c.gn.m = a.gn.m - b.gn.m;
image(a.gn.m);
image(b.gn.m);
image(c.gn.m);

## Conclusion: MI signals near or above the upper end of the noise level (here 0.2)
##   are detectable, but those near the median noise (here 0.04) disappear.
## It would be therefore useful to determine the median noise level and to
##   specify a signal-to-noise cut-off accordingly.

#_______________________________________________________________________________
## What is the distribution of MI values for a random matrix?
## compute MI values for given matrix
comp_MI = function(rmat, mimat) {
	cmb = combn(c(1:dim(rmat)[2]), 2);
	for (k in 1:dim(cmb)[2]) {
		i = cmb[1, k];
		j = cmb[2, k];
		x2d = rbind(rmat[ , i], rmat[, j]);
		mimat[i, j] = mi.empirical(x2d);
		mimat[j, i] = mimat[i, j];
	};
	return(mimat);
}

r = abs(rnorm(10000));
r = (r - min(r)) / (max(r) - min(r));
r.m = r;
dim(r.m) = c(100, 100);
mi.m = matrix(0, nrow = dim(r.m)[2], ncol = dim(r.m)[2]);

MI.m = comp_MI(r.m, mi.m);
MI.1.m = MI.m;
dim(MI.1.m) = c(dim(MI.m)[1] * dim(MI.m)[2], 1);

hist(MI.1.m[ ,1]);

mean(MI.1.m[,1]);
sd(MI.1.m[,1]);

## Observation: The median of the MI values is close to
##   the median of the random matrix.

## Conclusion: The MI values of a random value matrix are normal distributed
##   and sd decreases with increasing sample size (matrix dimension).
##   Therefore we can derive an error model from the MI matrix values.

#===============================================================================

