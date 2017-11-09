#! /usr/bin/R
#===============================================================================
#### Superimposition of the PCA loadings of Molecular Dynamics trajectories ####
### James Macpherson and Jens Kleinjung, 2016 (The Francis Crick Institute) ###
#===============================================================================

## Load the necessary R packages into your environment
library("bio3d");
setup.ncore(4, bigmem = FALSE);
library("ggplot2");
library('ggbiplot')
library('cluster')
library('NbClust')
library('ggrepel')
#______________________________________________________________________________
# Command line arguments
#______________________________________________________________________________

args = commandArgs(TRUE);
print('################################################################')
print('Usage: Rscript pca.R <reference.pdb> <projection.pdb> <output prefix>');
print('<reference.pdb>: PDB trajectory file used as the reference in PCA analysis.')
print('<projection.pdb>: PDB trajectory file projected into PCA space of the reference trajectory.')
print('<output prefix>: Output prefix for plots')

if(length(args) < 3)
{
	  stop('Not enough arguments. Please supply 3 arguments.')
}

#______________________________________________________________________________
# Read PDB files and adjust the trajectory dimensions
#______________________________________________________________________________

## Load reference and projected trajectories
## Change these variables to include the trajectories you wish to load!
traj_ref = args[1]   # path to the directory of the reference trajectory
traj_proj = args[2]  # path to the directory of the projected trajectory


ref = read.pdb(traj_ref, multi=TRUE)

proj = read.pdb(traj_proj, multi=TRUE)

## Fit conformers to remove roto-translations
xyz.ref = fit.xyz(fixed = ref$xyz, mobile = ref)
xyz.proj = fit.xyz(fixed = proj$xyz, mobile = proj)

# Change the dimensions of the trajectory matrix so that the reference and projected trajectories
# have the same number of residues.
dim.ref = dim(xyz.ref)[2]
dim.proj = dim(xyz.proj)[2]
diff.dim = dim.ref-dim.proj


if(diff.dim < 0){
    diff.dim = diff.dim * -1
}
if(dim.ref > dim.proj){
    xyz.ref = xyz.ref[ ,-c(1:diff.dim)]
}   	else if (dim.ref < dim.proj){
    xyz.proj = xyz.proj[ ,-c(1:diff.dim)]
}

#______________________________________________________________________________
# PCA analysis of the trajectories
#______________________________________________________________________________

# PCA plotting funtion
PCAggplot = function(pca.dat, pcx, pcy, ggplot.colours) {
    theme_set(theme_gray(base_size = 16));
    ggplot(data = as.data.frame(pca.dat$x), aes_string(x = pcx, y = pcy)) +
    geom_hline(yintercept = 0, colour = "gray65") +
    geom_vline(xintercept = 0, colour = "gray65") +
    geom_point(size = 4, colour = ggplot.colours);
}

# cumulative variance function
cumvar = function(pca.dat, ...){
	vars = apply(pca.dat$x, 2, var)
	props = vars / sum(vars)

	plot(cumsum(props),
	     type='o',
	     xlab='Principal Component',
	     ylab='Cumulative Variance (%)')
}

## Calculate the PC loadings for the reference trajectory
pca.ref = prcomp(xyz.ref, retx = TRUE, center = TRUE, scale. = TRUE);

## Calculate the PC loadings for the projected trajectory
pca.proj = prcomp(xyz.proj, retx = TRUE, center = TRUE, scale. = TRUE); #changed center to FALSE

## plot individual PCA analyses
pdf(paste(args[3], '_cumvar.pdf', sep=''))
par(mfrow=c(2,2))
cumvar(pca.ref, main='Reference trajectory')
cumvar(pca.proj, main='Projected trajectory')
dev.off()



## Project the projected trajectory PCA into the reference
projection.proj = predict(pca.ref, xyz.proj);



## Plot the projected and reference PCA on the same plot
pcascatterplt = function(pca.ref.dat.x, pca.ref.dat.y, pca.proj.dat.x, pca.proj.dat.y){

	# arrange the data frames for the PCA plot of the reference and projected
	# trajectories
	ref.dat = as.data.frame(cbind(pca.ref.dat.x, pca.ref.dat.y))
	names(ref.dat) = c('V1', 'V2')
	proj.dat = as.data.frame(cbind(pca.proj.dat.x, pca.proj.dat.y))
	names(proj.dat) = c('V1', 'V2')
 
	# define the dimensions of the projected PCA plot
	ref.dat.range = range(ref.dat)
	proj.dat.range = range(proj.dat)

	if(ref.dat.range[1] < proj.dat.range[1]){
		dim.min = ref.dat.range[1] - (ref.dat.range[1] * 0.05)
	}	else {
		dim.min = proj.dat.range[1] - (proj.dat.range[1] * 0.05)
	}

	if(ref.dat.range[2] > proj.dat.range[2]){
		dim.max = ref.dat.range[2] + (ref.dat.range[2] * 0.05)
	}	else {
		dim.max = proj.dat.range[2] + (proj.dat.range[2] * 0.05)
	}

	# plot the projection PCA
	plt1 = ggplot(data=ref.dat, aes(x=V1, y=V2)) +
     		geom_point(col='red', alpha=0.6, size=3) +
		xlim(dim.min, dim.max) +
		ylim(dim.min, dim.max) +
		xlab('PC1') +
		ylab('PC2') +

		geom_point(data=proj.dat, aes(x=V1, y=V2),
			   col='forestgreen', alpha=0.6, size=3) +
		theme_bw(base_size = 20) 

	# save the projected PCA plot to a pdf file
	pdf(paste(args[3], '_projected_pca.pdf', sep=''))
	par(mar=c(5,5,2,2))
	print(plt1)
	dev.off()

	# cluster the distribution in PC1 and PC2
	pdf(paste(args[3], '_pca_cluster.pdf', sep=''))
	dat = as.data.frame(rbind(ref.dat, proj.dat))
	nc = NbClust(dat, min.nc=2, max.nc=15, method="kmeans")

	plt2 = barplot(table(nc$Best.n[1,]),
		xlab="Numer of Clusters",
		ylab="Number of Criteria")

	fit = kmeans(dat, 2)
	plt3 = clusplot(dat, fit$cluster,
		 color=TRUE,
		 shade=TRUE)

	dev.off()
}

pcascatterplt(pca.ref$x[,1], pca.ref$x[,2], projection.proj[,1], projection.proj[,2])


## Compute the residue loadings onto PC1 and PC2
pc.loadings = function(ref_traj, proj_traj, ...){
	pc.ref = pca.xyz(ref_traj)
	pc.proj = pca.xyz(proj_traj)

	pdf(paste(args[3], '_pc_loadings_plt.pdf', sep=''))
	par(mfrow=c(2,1))
	plot(seq(from=13, to=529),
		   pc.ref$au[,1],
		   ylab='PC1 (Å)',
		   xlab='Residue Position',
		   typ='l',
		   col='red',
		   ylim=c(0,0.105))
	
	points(seq(from=13, to=529),
	       pc.proj$au[,1],
	       col='forestgreen',
	       typ='l')
	
	plot(seq(from=13, to=529),
		   pc.ref$au[,2],
		   ylab='PC2 (Å)',
		   xlab='Residue Position',
		   typ='l',
		   col='red',
		   ylim=c(0,0.105))
	
	points(seq(from=13, to=529),
	       pc.proj$au[,2],
	       col='forestgreen',
	       typ='l')

	dev.off()


	# difference in the PC loadings
	load.dif.1 = as.data.frame(cbind(
					 seq(from=13, to=529),
					 pc.ref$au[,1] - pc.proj$au[,1]))
	
	load.dif.2 = as.data.frame(cbind(
					 seq(from=13, to=529),
					 pc.ref$au[,2] - pc.proj$au[,2]))
	
	## identify residues significantly contributing to PC1 and PC2 differences Apo vs. FBP
	thresholdplt = function(datin){
		names(datin) = c('V1', 'V2')
		pfdr = mean(datin$V2) + (2*sd(datin$V2))
		nfdr = mean(datin$V2) - (2*sd(datin$V2))
		datin$threshold = as.factor(datin$V2 > pfdr | datin$V2 < nfdr)

		plt = ggplot(datin, aes(V1, V2)) +
			geom_point(alpha=0.6, aes(colour=threshold), size=3) +
			scale_color_manual(values = c("grey", "orange")) +
			scale_x_continuous(expand = c(0, 0)) +
			scale_y_continuous(expand = c(0, 0)) +
			ylim(min(datin$V2), max(datin$V2)) +
			theme_bw() +
			theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5), 'lines'),
					 axis.text = element_text(size=rel(2), colour="black"),
					 axis.title = element_text(size=rel(2), colour="black"),
					 legend.position="none") +
			xlab("Residue") +
			ylab('Loadings difference') +

			geom_text_repel(subset(datin, V2 > pfdr | V2 < nfdr),
					mapping = aes(label = V1),
					size = 3,
					color = 'black',
					box.padding = unit(0.5, "lines"),
					point.padding = unit(0.5, "lines"))

		return(plt)
	}

	pdf(paste(args[3], '_diff_loadings.pdf', sep=''), height=5, width=7)
	print(thresholdplt(load.dif.1))
	print(thresholdplt(load.dif.2))
	
	dev.off()
}



pc.loadings(xyz.ref, xyz.proj)
