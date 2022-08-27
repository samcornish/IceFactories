% script to collect SLP from each grid point in each model, winter means. 
clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/PSL/
addpath /home/ocean_personal_data/samc/CESM/variables/

load fnames_div.mat

lat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.PSL.192001-200512.nc','lat');
lon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.PSL.192001-200512.nc','lon');

for k = 1:32	% total number of ensembles that we will use
% 20th CENTURY
	SLP = ncread(strcat(fnames_div_20C{k},'.cam.h0.PSL.192001-200512.nc'),'PSL');

for n = 1:84 % no of years - 2 incomplete winter seasons
    Wslp = zeros(size(SLP,1),size(SLP,2),7);
    WslpE = zeros(144,size(SLP,2),7);
    % take winter months
    Wslp = SLP(:,:,12*(n-1)+10:12*(n-1)+16);
    WslpE = Wslp(145:end,:,:); % take Eastern hemisphere
    mean_slpE_20C(:,:,n,k) = nanmean(WslpE,3);
end
    fprintf('finished 20C ensemble # %d\n',k)

    clearvars SLP 
% RCP8.5 
	SLP = ncread(strcat(fnames_div_RCP85{k},'.cam.h0.PSL.200601-208012.nc'),'PSL');

for n = 1:73 % no of years - 2 incomplete winter seasons
    Wslp = zeros(size(SLP,1),size(SLP,2),7);
    WslpE = zeros(144,size(SLP,2),7);
    % take winter months
    Wslp = SLP(:,:,12*(n-1)+10:12*(n-1)+16);
    WslpE = Wslp(145:end,:,:); % take Eastern hemisphere
    mean_slpE_RCP85(:,:,n,k) = nanmean(WslpE,3);
end

fprintf('finished RCP8.5 ensemble # %d\n',k)
end
save /home/ocean_personal_data/samc/CESM/variables/CESM_SLP.mat mean_slpE_20C mean_slpE_RCP85 lat lon
