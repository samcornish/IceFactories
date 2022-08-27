% script to calculate and save the large arrays of winter div and I2O from enough models to do the scatter plot ice factory comparison
clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_personal_data/samc/CESM/variables/

load KLmask.mat
load fnames_div.mat
load fnames_fresh.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;


for k = 1:13	% total number of ensembles that we will use
% 20th CENTURY
	div = ncread(strcat(fnames_div_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
	I2O = ncread(strcat(fnames_fresh_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh'); 

msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
for i = 1:size(div,3)
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end

for n = 70:80 % no of years - 2 incomplete winter seasons
    Wdiv = zeros(size(div,1),size(div,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    Wdiv_store(:,:,size(Wdiv,3)*(n-70)+1:size(Wdiv,3)*(n-69),k) = Wdiv;
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    WI2O_store(:,:,size(WI2O,3)*(n-70)+1:size(WI2O,3)*(n-69),k) = WI2O;
    clearvars Wdiv WI2O
end
    fprintf('finished 20C ensemble # %d\n',k)
end
save /home/ocean_personal_data/samc/CESM/variables/CESM_div_I2O_KL_largearrays.mat Wdiv_store WI2O_store
