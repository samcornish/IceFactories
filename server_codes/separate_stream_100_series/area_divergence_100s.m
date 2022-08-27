
% script to calculate the area diverged (total and net) in the KL seas. On both freezing and non-freezing days. 
clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/aice/

load KLmask.mat
load fnames40.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
%cellarea_msk1xx = cellarea1xx.*mask1xx; % the members 1xx have a different grid
cellarea_msk = cellarea.*mask;	% masked cellarea
rep_cellarea = repmat(cellarea,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise
%rep_cellarea1xx = repmat(cellarea1xx,[1,1,7]); 
for k = 1:35	% total number of ensembles that we will use
% 20th CENTURY
if k == 1
	div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.185001-200512.nc'),'divu',[1,1,841],[291,48,1032]);
        I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.185001-200512.nc'),'fresh',[1,1,841],[291,48,1032]);
        SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.185001-200512.nc'),'aice',[1,1,841],[291,48,1032])/100;

else
	div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
	I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh'); 
	SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.192001-200512.nc'),'aice')/100;
end
msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
msk_SIC = zeros(size(SIC));
for i = 1:size(div,3)
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
msk_SIC(:,:,i) = SIC(:,:,i).*mask;
end

for n = 1:85 % no of years - 1 incomplete winter season
    Wdiv = NaN*zeros(size(div,1),size(div,2),7);
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
    WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    WSIC = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+16);
    fp = find(WSIC>0.15); % find cells with more than 15% sea ice
    % time integrated area divergence (m^2) values
    Wareadiv = Wdiv.*rep_cellarea*(30/100); % multiplying through by area and converting from %/day to area in a month
    fd = find(Wdiv>0); % find divergent points
    total_areadiv_20C(n,k) = nansum(Wareadiv(fd)); % find total area diverged (ie. only the positive part)
    fn = find(~isnan(Wdiv)); % find non NaN values
    net_areadiv_20C(n,k) = nansum(Wareadiv(fn)); % find net area diverged (allowing convergence to cancel divergence)
    monthly_net_areadiv_20C(n,:,k) = squeeze(nansum(nansum(Wdiv.*rep_cellarea*30/100))); % find the monthly structure of net area diverged  
    ffrz = find(WI2O<0); % find freezing points in time and space
    fdf = intersect(fd,ffrz);
    fpd = intersect(fp,fd);	% points that are both divergent and ice is present at more than 15%
    total_frz_areadiv_20C(n,k) = nansum(Wareadiv(fdf)); % total area diverged on freezing days
   % time mean divergence values (/s)
    meanposdiv_20C(n,k) = sum((Wdiv(fpd).*rep_cellarea(fpd))/sum(rep_cellarea(fpd))); % find the mean of the positive divergence, using the area-weighted array. On points where SIC is > 15% 
    meanposdiv_frz_20C(n,k) = sum((Wdiv(fdf).*rep_cellarea(fdf))/sum(rep_cellarea(fdf))); % find the mean of the positive divergence, on freezing points only, using the area-weighted array 
    meandiv_20C(n,k) = sum((Wdiv(fn).*rep_cellarea(fn))/sum(rep_cellarea(fn))); % find mean of divergence, on non-NaN cells (should be all ocean cells) - this should be equivalent to mean area export
    netdiv_over_ice_area_20C(n,k) = sum((Wdiv(fn).*WSIC(fn).*rep_cellarea(fn))/sum(WSIC(fn).*rep_cellarea(fn))); % this is the net div over the area of the ice
    netdiv_over_ice_extent_20C(n,k) = sum((Wdiv(fp).*rep_cellarea(fp))/sum(rep_cellarea(fp))); % this is the net div over the cells where ice is present at > 15 %
    posdiv_over_ice_extent_20C(n,k) = sum((Wdiv(fpd).*rep_cellarea(fpd))/sum(rep_cellarea(fp))); % this is mean positive div, counting only positive div but normalising over the entire ice area where SIC > 15%
    netdiv_frz_20C(n,k) = sum((Wdiv(ffrz).*rep_cellarea(ffrz))/sum(rep_cellarea(ffrz))); % find mean of divergence, on freezing cells  
    posdiv_frz_20C(n,k) = sum((Wdiv(fdf).*rep_cellarea(fdf))/sum(rep_cellarea(ffrz))); % find mean of positive divergence on freezing cells, normalised by total area of freezing cells  

    clearvars Wdiv WI2O WSIC
end
    fprintf('finished 20C ensemble # %d\n',k)
% now save the first part of the crossover year
n = 86;
  Wdiv = NaN*zeros(size(div,1),size(div,2),7);
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
    WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
    % take winter months
    Wdiv(:,:,1:3) = msk_div(:,:,12*(n-1)+10:12*(n-1)+12);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O(:,:,1:3) = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    WSIC(:,:,1:3) = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+12);
    clearvars msk_I2O msk_div I2O div 

% RCP8.5 
if k >= 34
        div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-210012.nc'),'divu',[1,1,1],[291,48,900]);
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-210012.nc'),'fresh',[1,1,1],[291,48,900]);
        SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-210012.nc'),'aice',[1,1,1],[291,48,900])/100;
else
	div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-208012.nc'),'divu');
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-208012.nc'),'fresh');
	SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-208012.nc'),'aice')/100;
end
msk_SIC = zeros(size(SIC));
msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
for i = 1:size(div,3)
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
msk_SIC(:,:,i) = SIC(:,:,i).*mask;
end

% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
    Wdiv(:,:,4:7) = msk_div(:,:,1:4);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O(:,:,4:7) = msk_I2O(:,:,1:4).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    WSIC(:,:,4:7) = msk_SIC(:,:,1:4);
    fp = find(WSIC>0.15); % find cells with more than 15% sea ice
    Wareadiv = Wdiv.*rep_cellarea*(30/100); % multiplying through by area and converting from %/day to area in a month
    fd = find(Wdiv>0); % find divergent points
    total_areadiv_20C(n,k) = nansum(Wareadiv(fd)); % find total area diverged (ie. only the positive part)
    fn = find(~isnan(Wdiv)); % find non NaN values
    net_areadiv_20C(n,k) = nansum(Wareadiv(fn)); % find net area diverged (allowing convergence to cancel divergence)
    monthly_net_areadiv_20C(n,:,k) = squeeze(nansum(nansum(Wdiv.*rep_cellarea*30/100))); % find the monthly structure of net area diverged  
    ffrz = find(WI2O<0); % find freezing points in time and space
    fdf = intersect(fd,ffrz);
    fpd = intersect(fp,fd);     % points that are both divergent and ice is present at more than 15%
    total_frz_areadiv_20C(n,k) = nansum(Wareadiv(fdf)); % total area diverged on freezing days
   % time mean divergence values (/s)
    meanposdiv_20C(n,k) = sum((Wdiv(fpd).*rep_cellarea(fpd))/sum(rep_cellarea(fpd))); % find the mean of the positive divergence, using the area-weighted array. On points where SIC is > 15% 
    meanposdiv_frz_20C(n,k) = sum((Wdiv(fdf).*rep_cellarea(fdf))/sum(rep_cellarea(fdf))); % find the mean of the positive divergence, on freezing points only, using the area-weighted array 
    meandiv_20C(n,k) = sum((Wdiv(fn).*rep_cellarea(fn))/sum(rep_cellarea(fn))); % find mean of divergence, on non-NaN cells (should be all ocean cells) - this should be equivalent to mean area export
    netdiv_over_ice_area_20C(n,k) = sum((Wdiv(fn).*WSIC(fn).*rep_cellarea(fn))/sum(WSIC(fn).*rep_cellarea(fn))); % this is the net div over the area of the ice
    netdiv_over_ice_extent_20C(n,k) = sum((Wdiv(fp).*rep_cellarea(fp))/sum(rep_cellarea(fp))); % this is the net div over the cells where ice is present at > 15 %
    posdiv_over_ice_extent_20C(n,k) = sum((Wdiv(fpd).*rep_cellarea(fpd))/sum(rep_cellarea(fp))); % this is mean positive div, counting only positive div but normalising over the entire ice area where SIC > 15%
    netdiv_frz_20C(n,k) = sum((Wdiv(ffrz).*rep_cellarea(ffrz))/sum(rep_cellarea(ffrz))); % find mean of divergence, on freezing cells  
    posdiv_frz_20C(n,k) = sum((Wdiv(fdf).*rep_cellarea(fdf))/sum(rep_cellarea(ffrz))); % find mean of positive divergence on freezing cells, normalised by total area of freezing cells  


    clearvars Wdiv WI2O WSIC



for n = 1:73 % no of years - 2 incomplete winter seasons
    Wdiv = zeros(size(div,1),size(div,2),7);
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
    WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    WSIC = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+16);
    fp = find(WSIC>0.15); % find cells with more than 15% sea ice
    Wareadiv = Wdiv.*rep_cellarea*(30/100); % multiplying through by area and converting from %/day to area in a month
    fd = find(Wdiv>0); % find divergent points
    total_areadiv_RCP85(n,k) = nansum(Wareadiv(fd)); % find total area diverged (ie. only the positive part)
    fn = find(~isnan(Wdiv)); % find non NaN values
    net_areadiv_RCP85(n,k) = nansum(Wareadiv(fn)); % find net area diverged (allowing convergence to cancel divergence)
    monthly_net_areadiv_RCP85(n,:,k) = squeeze(nansum(nansum(Wdiv.*rep_cellarea*30/100))); % find the monthly structure of net area diverged  
    ffrz = find(WI2O<0); % find freezing points in time and space
    fdf = intersect(fd,ffrz);
    fpd = intersect(fp,fd);	% points that are both divergent and ice is present at more than 15%
    total_frz_areadiv_RCP85(n,k) = nansum(Wareadiv(fdf)); % total area diverged on freezing days
    % time mean divergence values (/s)
    meanposdiv_RCP85(n,k) = sum((Wdiv(fpd).*rep_cellarea(fpd))/sum(rep_cellarea(fpd))); % find the mean of the positive divergence, using the area-weighted array 
    meanposdiv_frz_RCP85(n,k) = sum((Wdiv(fdf).*rep_cellarea(fdf))/sum(rep_cellarea(fdf))); % find the mean of the positive divergence, on freezing points only, using the area-weighted array
    meandiv_RCP85(n,k) = sum((Wdiv(fn).*rep_cellarea(fn))/sum(rep_cellarea(fn))); % find mean of divergence, on non-NaN cells (should be all ocean cells) - this should be equivalent to mean area export
    netdiv_over_ice_area_RCP85(n,k) = sum((Wdiv(fn).*WSIC(fn).*rep_cellarea(fn))/sum(WSIC(fn).*rep_cellarea(fn))); % this is the net div over the area of the ice
    netdiv_over_ice_extent_RCP85(n,k) = sum((Wdiv(fp).*rep_cellarea(fp))/sum(rep_cellarea(fp))); % this is the net div over the cells where ice is > 15 % conc
    posdiv_over_ice_extent_RCP85(n,k) = sum((Wdiv(fpd).*rep_cellarea(fpd))/sum(rep_cellarea(fp))); % this is mean positive div, counting only positive div but normalising over the entire ice area where SIC > 15%
    netdiv_frz_RCP85(n,k) = sum((Wdiv(ffrz).*rep_cellarea(ffrz))/sum(rep_cellarea(ffrz))); % find mean of divergence, on freezing cells  
    posdiv_frz_RCP85(n,k) = sum((Wdiv(fdf).*rep_cellarea(fdf))/sum(rep_cellarea(ffrz))); % find mean of positive divergence on freezing cells, normalised by total area of freezing cells  

    clearvars Wdiv WI2O WSIC
end
    fprintf('finished RCP85 ensemble # %d\n',k)
end



save /home/ocean_personal_data/samc/CESM/variables/CESM_areadiv_20C_RCP85.mat total_areadiv_20C net_areadiv_20C monthly_net_areadiv_20C total_frz_areadiv_20C... 
total_areadiv_RCP85 net_areadiv_RCP85 monthly_net_areadiv_RCP85 total_frz_areadiv_RCP85 meanposdiv_20C meanposdiv_RCP85 meanposdiv_frz_20C meanposdiv_frz_RCP85 meandiv_20C meandiv_RCP85 netdiv_over_ice_area_20C netdiv_over_ice_area_RCP85 netdiv_over_ice_extent_20C netdiv_over_ice_extent_RCP85 netdiv_frz_20C netdiv_frz_RCP85 posdiv_over_ice_extent_20C posdiv_over_ice_extent_RCP85 posdiv_frz_20C posdiv_frz_RCP85
