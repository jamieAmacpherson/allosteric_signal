#! /usr/bin/R

dat = read.table("all_rms.dat", header = TRUE);

dat.apo1 = dat[ , "bound1"] == "apo";
dat.apo2 = dat[ , "bound2"] == "apo";
dat.fbp1 = dat[ , "bound1"] == "fbp";
dat.fbp2 = dat[ , "bound2"] == "fbp";

dat.apo1.apo2 = dat.apo1 & dat.apo2;
hist(dat[dat.apo1.apo2, "RMSD"]);

dat.fbp1.fbp2 = dat.fbp1 & dat.fbp2;
hist(dat[dat.fbp1.fbp2, "RMSD"]);

dat.apo1.fbp2 = dat.apo1 & dat.fbp2;
dat.fbp1.apo2 = dat.fbp1 & dat.apo2;
dat.mixed = dat.apo1.fbp2 | dat.fbp1.apo2;
hist(dat[dat.mixed, "RMSD"]);



