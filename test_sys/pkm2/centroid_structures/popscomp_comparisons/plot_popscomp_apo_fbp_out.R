#! /usr/bin/R

#===============================================================================
# plot POPSCOMP output for apo and fbp
#===============================================================================

## output 
inFileName = list.files(path = ".", full.names = TRUE, pattern = 'out$');

data.l = lapply(inFileName, function(x){read.table(file = x, header = TRUE)});
data.df = do.call(rbind.data.frame, data.l);

data.df = cbind(data.df, rep(substring(inFileName, 3, 7), each = 6));
colnames(data.df) = c("domainPair", "APO.Phob", "APO.Phil", "APO.Total",
									"FBP.Phob", "FBP.Phil", "FBP.Total",
					  "source");

## plot data
## AB
pdf("deltaSASA_AB.pdf");
boxplot(data.df[data.df$"domainPair" == "Adom.Bdom", 2:7],
		col = c("green", "blue", "black", "green", "blue", "black"),
		ylab = "delta SASA / A^2");
dev.off();

## AC
pdf("deltaSASA_AC.pdf");
boxplot(data.df[data.df$"domainPair" == "Adom.Cdom", 2:7],
		col = c("green", "blue", "black", "green", "blue", "black"),
		ylab = "delta SASA / A^2");
dev.off();

## BC 
## practically no differences

## NA
pdf("deltaSASA_NA.pdf");
boxplot(data.df[data.df$"domainPair" == "Ndom.Adom", 2:7],
		col = c("green", "blue", "black", "green", "blue", "black"),
		ylab = "delta SASA / A^2");
dev.off();

## NB
## practically no differences

## NC 
pdf("deltaSASA_NC.pdf");
boxplot(data.df[data.df$"domainPair" == "Ndom.Cdom", 2:7],
		col = c("green", "blue", "black", "green", "blue", "black"),
		ylab = "delta SASA / A^2");
dev.off();


