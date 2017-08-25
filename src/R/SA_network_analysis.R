######################################################################
#
#       process data
#
######################################################################
source('/home/macphej/jm.software/apps/gsatools-4.5.x-1.00/Rscripts/modules/FragMINetwork_JM_ver_20160527.R')

######################################################################
# Imports
library(ggplot2)
library(ggrepel)
library(optparse)

######################################################################
#       general parameters 
fragLength = 4
FDRcutoff = 0.001
percCutoff = 0.0
sigmas = 2 		# number of standard deviations as signficance
			# cutoff
######################################################################
#       graphical parameters 
vertexColor = "#55555555"
vertexFrameColor = "#000000ff"
edgeColor = "#55555522"

altVertexLabelColor = "#ffffffff"
altVertexColor = "#0000ffff"
altVertexFrameColor = "#0000ffff"
altEdgeColor = "#ff00ff77"

nColor = 50

######################################################################
#
# Parameters
#
# The script can easily be adapted to process user-provided data
# by simply modifying the values of these parameters.
#
wSize = 1200
hSize = 800

# command line flags
option_list = list(
	make_option(c("-pdb", "--PDB_file"), type="character", default=NULL, 
              help="pdb file name", metavar="character"),

	make_option(c("-mi", "--mutual_information"), type="character", default="NULL", 
              help="Mutual information matrix file name", metavar="character"),

	make_option(c("-je", "--joint_entropy"), type="character", default="NULL", 
              help="Joint entropy matrix file name", metavar="character"),

	make_option(c("-err", "--error_est"), type="character", default="NULL", 
              help="Estimated error matrix file name", metavar="character"),

	make_option(c("-pval", "--p_value"), type="character", default="NULL", 
              help="P value matrix file name", metavar="character"),
	
	make_option(c("-o", "--output"), type="character", default="NULL", 
              help="Output prefix", metavar="character")
	
); 

# parse command line arguments 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# exit if PDB file is not supplied
if (is.null(opt$pdb)){
  print_help(opt_parser)
  stop("Arguments must be supplied.", call.=FALSE)
}

# assign file names
pdbFilename = opt$PDB_file
MIInputFilename = opt$mutual_information
jointEntropyInputFilename = opt$joint_entropy
estErrorInputFilename = opt$error_est
pValueInputFilename = opt$p_value

outputPrefix = opt$output

######################################################################
#
# load data
loadmat = function(filename){
	return(as.matrix(read.table(filename)))
}

MILF = loadmat(MIInputFilename)
jHLF = loadmat(jointEntropyInputFilename)
eMILF = loadmat(estErrorInputFilename)
pMILF = loadmat(pValueInputFilename)

pdbDm = calculateDistanceMatrix(pdbFilename, fragLength) 

nFrag = dim(MILF)[1]


######################################################################
#
# Network graph
#
######################################################################

enpdLFMatrix = filterMatrixContactLongrange(pdbDm,
			filterMIMatrix(renameMatrix(MILF),
					pMILF,
					eMILF,
					jHLF,
					FDRcutoff),
					TRUE)

write.table(enpdLFMatrix,
			file = paste(outputPrefix, '_enpdMI.mat', sep = ''),
					quote = F,
					row.names =F,
					col.names = F)

write.table(pdbDm,
		file = paste(outputPrefix, '_dismat.mat', sep = ''),
			quote = F,
			row.names = F,
			col.names = F)


LFgs = returnGraph(pdbFilename, fragLength, 
                     MILF, pMILF, eMILF, jHLF,
                     FDRcutoff)

pcalayout = layout.norm(princomp(pdbDm)$scores[,c(1,2)], -1,1, -1,1)

outpcalayout = pcalayout
row.names(outpcalayout) = V(LFgs)$name
write.table(outpcalayout, file = paste(outputPrefix, '.pcalayout', sep = ''), quote = F, col.names = F)

write.graph(LFgs, file = paste(outputPrefix, '.graphml', sep = ''), format= "graphml")
write.graph(LFgs, file = paste(outputPrefix, '.gml', sep = ''), format= "gml")

reWeightVector = 1 - (E(LFgs)$weight - min(E(LFgs)$weight)) / (max(E(LFgs)$weight) - min(E(LFgs)$weight))
weightColorVector = rgb.palette(nColor)[as.integer(reWeightVector * nColor) + 1]
weightColorData = data.frame(get.edgelist(LFgs, names=TRUE), wc = weightColorVector)
write.table(weightColorData, file = paste(outputPrefix, '.colorData', sep = ''), quote = F, row.names = F, col.names = F)

system(paste('python /home/macphej/jm.software/apps/gsatools-4.5.x-1.00/Rscripts/modules/igraph2cytoscape.py ',
	outputPrefix,
	'.pcalayout ',
	outputPrefix,
	'.gml ',
	outputPrefix,
	'.colorData',
	sep =''))


######################################################################
# Mutual information distribution plot function
######################################################################

# Normalisation function
normalisedat = function(dat){
	tmp = (dat - min(dat)) / (max(dat) - min(dat))
	return(tmp)
} 

# plot the mutual information dist function
MIdistribution = function(nmi_matrix, difference){
	
	# set infinite entries in the matrix to zero
	nmi_matrix[is.infinite(nmi_matrix)] <- 0
	
	if(max(nmi_matrix)>1){
		nmi_matrix = normalisedat(nmi_matrix)
	}
	
	mathist = hist(nmi_matrix, plot=FALSE, breaks=2000)
	pfdr = mean(nmi_matrix) + (sigmas*sd(nmi_matrix))
	nfdr = mean(nmi_matrix) - (sigmas*sd(nmi_matrix))

	plot(mathist$breaks[-1],
		log(mathist$counts),
		type='h',
		xlab = expression(italic(nMI)),
		ylab = expression(ln ~ italic(f(nMI))),
		cex.lab = 2,
		cex.axis = 2,
		panel.first = grid(),
		ylim = c(0, 16))

}

pdf(paste(outputPrefix, "_nMI_dist.pdf", sep=""))
par(mar=c(5,5,1,1))
par(mfrow=c(3,1))
MIdistribution(MIInputFilename, FALSE)
MIdistribution(enpdLFMatrix, FALSE)
dev.off()

######################################################################
# Eigenvector centrality plotting function
######################################################################
plotevcent = function(evcentdat){
	evcent = as.numeric(unlist(evcentdat))
        mathist = hist(evcent, plot=FALSE, breaks=2000)
        pfdr = mean(evcent) + (4*sd(evcent))

	pdf("evdist.pdf")
	par(mar=c(5,5,1,1))
        plot(mathist$breaks[-1],
                mathist$counts,
                type='h',
        #       log='y',
                xlab = expression(italic(x[nu])),
                ylab = expression(italic(f(x[nu]))),
                cex.lab = 2,
                cex.axis = 2,
                panel.first = grid())

        # line at mu + 3sd      
        abline(v=pfdr,
                lty=2,
                col="red",
                lwd=3)

        # text
        text(pfdr + 0.1, 30, expression(2*sigma), cex=2)
	dev.off()

	dat = as.data.frame(cbind(
			seq(from=1, to=nrow(evcentdat)),
			evcentdat))

	names(dat) = c("V1", "V2")

	dat$threshold = as.factor(dat$V2 > pfdr)
	
	pl = ggplot(dat, aes(V1, V2)) +
                geom_point(alpha=0.6, aes(colour=threshold)) +
                scale_color_manual(values = c("grey", "forestgreen")) +
                scale_x_continuous(expand = c(0, 0)) +
                scale_y_continuous(expand = c(0, 0)) +
                ylim(min(dat$V2), max(dat$V2)) +
                theme_bw() +
                theme(plot.margin = unit(c(0.2, 0.2, 0.5, 0.5),
                                "lines"),
                        axis.text = element_text(size=rel(2),
                                colour="black"),
                        axis.title = element_text(size=rel(2),
                                colour="black"),
                legend.position="none") +
                xlab("Fragment") +
                ylab(expression(italic(x[nu]))) +

            geom_text_repel(subset(dat, V2 > pfdr),
                       mapping = aes(label = V1),
                       size = 3,
                       color = 'black',
                       box.padding = unit(0.5, "lines"),
                       point.padding = unit(0.5, "lines"))
	
	pdf(paste(outputPrefix, "ev_plot.pdf", sep=""), height=5, width=7)	
	print(pl)	
	dev.off()

	write.table(evcentdat,
		file = paste(outputPrefix, '.evcent.dat' , sep = ''),
		row.names = F,
		col.names = F,
		quote = F)
	
        significant_frags = subset(dat, V2 > pfdr)$V1
        write.table(file=paste("sig_frags.txt", sep=""), t(significant_frags),
                row.names=FALSE,
                col.names=FALSE)
	
}
evcentdat = evcent(LFgs)$vector

plotevcent(evcentdat)
#


