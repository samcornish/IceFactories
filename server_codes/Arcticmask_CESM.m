% mask of central Arctic for CESM

clear
close all
addpath /home/ocean_shared_data2/CESM_arctic/divu/
load coastlines

lat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','TLAT');
lon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','TLON');

klat = find(lat<80);
klonE = find(lon<60);
klonW = find(lon>240);

kmaskE = intersect(klat,klonE);
kmaskW = intersect(klat,klonW);

ArcticMask = ones(291,48);
ArcticMask(kmaskE)=NaN;
ArcticMask(kmaskW)=NaN;



save /home/ocean_personal_data/samc/CESM/variables/ArcticMask.mat ArcticMask lat lon


