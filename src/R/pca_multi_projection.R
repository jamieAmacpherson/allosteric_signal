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
library('plot3D')
library('RColorBrewer')
library('plotly')
library('webshot')
library('scatterplot3d')
#______________________________________________________________________________
# Command line arguments
#______________________________________________________________________________

args = commandArgs(TRUE);
print('################################################################')
print('Usage: Rscript pca.R <reference.pdb> <projection1.pdb> <projection2.pdb> <projection3.pdb> <projection4.pdb> <output prefix>');
print('<reference.pdb>: PDB trajectory file used as the reference in PCA analysis.')
print('<projection1.pdb>: PDB trajectory file projected into PCA space of the reference trajectory.')
print('<projection2.pdb>: PDB trajectory file projected into PCA space of the reference trajectory.')
print('<projection3.pdb>: PDB trajectory file projected into PCA space of the reference trajectory.')
print('<projection4.pdb>: PDB trajectory file projected into PCA space of the reference trajectory.')
print('<number of c-alpha>: Number of C-alpha atoms to be projected (must be consistent across all trajectories)')
print('<output prefix>: Output prefix for plots')

if(length(args) < 7)
{
	  stop('Not enough arguments. See printed help.')
}

#______________________________________________________________________________
# Read PDB files and adjust the trajectory dimensions
#______________________________________________________________________________

load_traj = function(traj, atomdim){

	# read in the pdb file
	traj.pdb = read.pdb(traj, multi=TRUE)

	# remove rot-trans motions
	xyz.coords = fit.xyz(fixed = traj.pdb$xyz, mobile=traj.pdb)

	# determine the dimensions of the trajectory file
	targetdim = as.numeric(atomdim) * 3
	
	dim.traj = dim(xyz.coords)[2]

	# if the dimensions of the trajectory are greater than
	# the pre-defined limits, trim down the trajectory	
#	stopifnot(dim.traj != targetdim)

	return(xyz.coords)
}

## Load reference and projected trajectories
xyz.ref   = load_traj(args[1], args[6])  
xyz.proj1 = load_traj(args[2], args[6])  
xyz.proj2 = load_traj(args[3], args[6])  
xyz.proj3 = load_traj(args[4], args[6])  
xyz.proj4 = load_traj(args[5], args[6])  


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

# eigenanalysis of individual trajectories
traj.pca = function(trajectory_coordinates, ...){
	pca.dat = prcomp(trajectory_coordinates,
			 retx = TRUE,
			 center = TRUE,
			 scale. = TRUE)

	return(pca.dat)
}

## plot individual PCA analyses
pdf(paste(args[3], '_cumvar.pdf', sep=''))
par(mfrow=c(2,3))
cumvar(traj.pca(xyz.ref), main='Reference trajectory')
cumvar(traj.pca(xyz.proj1), main='Projected trajectory 1')
cumvar(traj.pca(xyz.proj2), main='Projected trajectory 2')
cumvar(traj.pca(xyz.proj3), main='Projected trajectory 3')
cumvar(traj.pca(xyz.proj4), main='Projected trajectory 4')
dev.off()



## Project the projected trajectory PCA into the reference
pca.projection = function(reference.coords, projection.coords){
	# pca of reference trajectory
	pca.ref = traj.pca(reference.coords)
	
	# apply the same rotation matrix to the projected 
	# trajectory
	projected.pca = predict(pca.ref, projection.coords)

	return(projected.pca)
}


## Plot the projected and reference PCA on the same plot
pcascatterplt = function(pca.ref.dat.x, pca.ref.dat.y, pca.ref.dat.z, pca.proj1.dat.x, pca.proj1.dat.y, pca.proj1.dat.z, pca.proj2.dat.x, pca.proj2.dat.y, pca.proj2.dat.z, pca.proj3.dat.x, pca.proj3.dat.y, pca.proj3.dat.z, pca.proj4.dat.x, pca.proj4.dat.y, pca.proj4.dat.z){

	# arrange the data frames for the PCA plot of the reference and projected
	# trajectories
	arrange.dat = function(xdat, ydat, zdat){
		tmp = as.data.frame(cbind(xdat, ydat, zdat))
		names(tmp) = c('V1', 'V2', 'V3')

		return(tmp)
	}

	ref.dat = arrange.dat(pca.ref.dat.x, pca.ref.dat.y, pca.ref.dat.z)
	
	proj1.dat = arrange.dat(pca.proj1.dat.x, pca.proj1.dat.y, pca.proj1.dat.z)
	proj2.dat = arrange.dat(pca.proj2.dat.x, pca.proj2.dat.y, pca.proj2.dat.z)
	proj3.dat = arrange.dat(pca.proj3.dat.x, pca.proj3.dat.y, pca.proj3.dat.z)
	proj4.dat = arrange.dat(pca.proj4.dat.x, pca.proj4.dat.y, pca.proj4.dat.z)
 
	# define the dimensions of the projected PCA plot
	list.dat = list(ref.dat, proj1.dat, proj2.dat, proj3.dat, proj4.dat)

	dim.min = range(list.dat)[1]
	dim.max = range(list.dat)[2]

	mycols <- c("tFBP", "tApo", "mFBP", "mApo", "mTEPP")

	# plot the projection PCA
	plt1 = ggplot(data=ref.dat, aes(x=V1, y=V2)) +
     		geom_point(aes(col=mycols[1]), alpha=0.6, size=3) +
		xlim(dim.min, dim.max) +
		ylim(dim.min, dim.max) +
		xlab('PC1') +
		ylab('PC2') +

		geom_point(data=proj1.dat, aes(x=V1, y=V2, col=mycols[2]),
			   alpha=0.6, size=3) +

		geom_point(data=proj2.dat, aes(x=V1, y=V2, col=mycols[3]),
			   alpha=0.6, size=3) +

		geom_point(data=proj3.dat, aes(x=V1, y=V2, col=mycols[4]),
			   alpha=0.6, size=3) +

		geom_point(data=proj4.dat, aes(x=V1, y=V2, col=mycols[5]),
			   alpha=0.6, size=3) +

		theme_bw(base_size = 20) +
	        
		scale_fill_manual(
				  name="Simulation",
				  labels=c("tFBP", "tApo", "mFBP", "mApo", "mTEPP"))		

	# save the projected PCA plot to a pdf file
	pdf(paste(args[7], '_projected_pca.pdf', sep=''))
	par(mar=c(5,5,2,2))
	print(plt1)
	dev.off()

	## 3d scatter-plot
	# function for adding points to 3d scatter plot
	sysstring = function(system, len){
		# initiate a data frame of strings 
		tmp = as.data.frame(rep(system, len))
		# give the data frame a consistent name so that
		# binding it by row doesnt throw up errors
		names(tmp) = 'sys'
		return(tmp)
	}

	# bind together all of the PCA projections and include a system name
	# so that the colours can be assigned automatically
	dat = as.data.frame(rbind(cbind(ref.dat, sysstring('tFBP', nrow(ref.dat))),
				  cbind(proj1.dat, sysstring('tApo', nrow(proj1.dat))),
				  cbind(proj2.dat, sysstring('mFBP', nrow(proj2.dat))),
				  cbind(proj3.dat, sysstring('mApo', nrow(proj3.dat))),
				  cbind(proj4.dat, sysstring('mTEPP', nrow(proj4.dat)))))
	names(dat) = c('x', 'y', 'z', 'sys')

	p = plot_ly(dat, x = ~x, y = ~y, z = ~z,
		    color = ~sys) %>%
			    layout(scene = list(xaxis = list(title = 'PC1'),
						yaxis = list(title = 'PC2'),
						zaxis = list(title = 'PC3'))) %>%
							add_markers()
	
	# save the plot to an HTML file
	htmlwidgets::saveWidget(p, "pca_plot_3d.html")

	# save the 3d plot to a png exported locally
	Sys.setenv("plotly_username" = "macpherson.js93")
	Sys.setenv("plotly_api_key" = "27fd4Bl1DFDj33Adrkiq")

	plotly_IMAGE(p, format = "png", out_file = 'pca_plot_3d.png')

	colors <- brewer.pal(n = 5, name = "Dark2")
	colors <- colors[as.numeric(dat$sys)]

	pdf('pca_plot_3d.pdf')
	scatterplot3d(dat$x,
		      dat$y,
		      dat$z,
		      color=colors,
		      xlab = 'PC1',
		      ylab = 'PC2',
		      zlab = 'PC3',
		      pch = 16)

	#legend("bottom", legend = levels(dat$sys),
	#             col =  colours, pch = 16, 
	#	     inset = -0.25, xpd = TRUE, horiz = TRUE)
	dev.off()

	# cluster the distribution in PC1 and PC2
	pdf(paste(args[3], '_pca_cluster.pdf', sep=''))
	dat = as.data.frame(rbind(ref.dat,
				  proj1.dat,
				  proj2.dat,
				  proj3.dat,
				  proj4.dat))

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

pcascatterplt(traj.pca(xyz.ref)$x[,1],
	      traj.pca(xyz.ref)$x[,2], 
	      traj.pca(xyz.ref)$x[,3],

	      pca.projection(xyz.ref, xyz.proj1)[,1],
	      pca.projection(xyz.ref, xyz.proj1)[,2],
	      pca.projection(xyz.ref, xyz.proj1)[,3],

	      pca.projection(xyz.ref, xyz.proj2)[,1],
	      pca.projection(xyz.ref, xyz.proj2)[,2],
	      pca.projection(xyz.ref, xyz.proj2)[,3],

	      pca.projection(xyz.ref, xyz.proj3)[,1],
	      pca.projection(xyz.ref, xyz.proj3)[,2],
	      pca.projection(xyz.ref, xyz.proj3)[,3],

	      pca.projection(xyz.ref, xyz.proj4)[,1],
	      pca.projection(xyz.ref, xyz.proj4)[,2],
	      pca.projection(xyz.ref, xyz.proj4)[,3])

## Compute the residue loadings onto PC1 and PC2

pc.loadings = function(ref_traj, ref_name, proj_traj, proj_name, outputprefix, ...){
	pc.ref = pca.xyz(ref_traj)
	pc.proj = pca.xyz(proj_traj)

	# determine the range of the two trajectories for the loadings
        # along the first eigenvector	
	dat.list1 = list(pc.ref$au[,1], pc.proj$au[,1])
	dim1.min = range(dat.list1)[1]
	dim1.max = range(dat.list1)[2]
	
	# determine the range of the two trajectories for the loadings
        # along the second eigenvector	
	dat.list2 = list(pc.ref$au[,2], pc.proj$au[,2])
	dim2.min = range(dat.list2)[1]
	dim2.max = range(dat.list2)[2]

	#____________________________________________________
	## generate loadings plot along the first eigenvector
	#____________________________________________________
	pdf(paste(outputprefix, '_pc_loadings_plt.pdf', sep=''))
	par(mar=c(5,5,2,2))
	par(mfrow=c(2,1))
	plot(seq(from=13, to=530, length.out=length(pc.ref$au[,1])),
		   pc.ref$au[,1],
		   ylab='PC1 (Å)',
		   xlab='Residue Position',
		   typ='l',
#		   col='red',
		   ylim=c(dim1.min, dim1.max),
		   cex.axis=2,
		   cex.lab=2,
		   lty=1)
	
	points(seq(from=13, to=530, length.out=length(pc.proj$au[,1])),
	       pc.proj$au[,1],
#	       col='forestgreen',
	       typ='l',
	       lty=4)

	legend('topleft',
	       c(ref_name, proj_name),
	       lty = c(1,4))
	       #col = c('forestgreen', 'red'))
	#____________________________________________________
	## generate loadings plot along the second eigenvector
	#____________________________________________________
	plot(seq(from=13, to=530, length.out=length(pc.ref$au[,2])),
		   pc.ref$au[,2],
		   ylab='PC2 (Å)',
		   xlab='Residue Position',
		   typ='l',
		   lty=1,
		 #  col='red',
		   ylim=c(dim2.min, dim2.max),
		   cex.axis = 2,
		   cex.lab = 2)
	
	points(seq(from=13, to=530, length.out=length(pc.proj$au[,2])),
	       pc.proj$au[,2],
	       #col='forestgreen',
	       typ='l',
	       lty=4)
	
	legend('topleft',
	       c(ref_name, proj_name),
	       lty = c(1,4))
	       #col = c('forestgreen', 'red'))

	dev.off()


	# difference in the PC loadings
	load.dif.1 = as.data.frame(cbind(
					 seq(from=13, to=529),
					 pc.ref$au[,1] - pc.proj$au[,1]))
	
	load.dif.2 = as.data.frame(cbind(
					 seq(from=13, to=529),
					 pc.ref$au[,2] - pc.proj$au[,2]))
	
	## identify residues significantly contributing to PC1 and PC2 differences Apo vs. FBP
	thresholdplt = function(datin, PCn){
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
			xlab('Residue') +
			ylab(paste(PCn, ' loadings difference', sep='')) +

			geom_text_repel(subset(datin, V2 > pfdr | V2 < nfdr),
					mapping = aes(label = V1),
					size = 3,
					color = 'black',
					box.padding = unit(0.5, "lines"),
					point.padding = unit(0.5, "lines"))

		return(plt)
	}

	pdf(paste(outputprefix, '_diff_loadings.pdf', sep=''), height=5, width=7)
	print(thresholdplt(load.dif.1, 'PC1'))
	print(thresholdplt(load.dif.2, 'PC2'))
	
	dev.off()
}


#pc.loadings(xyz.ref, 'tFBP' , xyz.proj1, 'tApo', 'tFBP-tApo')
#pc.loadings(xyz.ref, 'tFBP', xyz.proj2, 'mFBP', 'tFBP-mFBP')
#pc.loadings(xyz.proj1, 'tApo', xyz.proj3, 'mApo', 'tApo-mApo')
#pc.loadings(xyz.proj2, 'mFBP', xyz.proj4, 'mTepp-46', 'mFBP-mTEPP')


