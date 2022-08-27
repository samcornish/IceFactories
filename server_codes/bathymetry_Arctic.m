% script to record the ocean bathymetry in the Arctic Ocean in CESM

clear
close all
addpath /home/ocean_shared_data2/CESM_arctic/TEMP/

% variabe HT is ocean depth at T points, in cm

bathymetry = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.TEMP.192001-200512.nc','HT')/100;
Tlat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.TEMP.192001-200512.nc','TLAT');
Tlon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.TEMP.192001-200512.nc','TLONG');

save /home/ocean_personal_data/samc/CESM/variables/bathymetry_Arctic_ocean.mat bathymetry Tlat Tlon


