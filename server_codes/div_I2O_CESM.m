% script to find the regression coefficients, variance explained, mean(+div), std(div), sum of I2O

clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_personal_data/samc/CESM/variables/

load KLmask.mat
load fnames40.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
cellarea_msk = cellarea.*mask;	% masked cellarea
rep_cellarea = repmat(cellarea,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise

for k = 1:35	% total number of ensembles that we will use
% 20th CENTURY
if k == 1
        div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.185001-200512.nc'),'divu',[1,1,841],[291,48,1032]);
        I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.185001-200512.nc'),'fresh',[1,1,841],[291,48,1032]);
else
	div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
	I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh'); 
end
msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
for i = 1:size(div,3)
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end

for n = 1:85 % no of years - 1 incomplete winter season
    Wdiv = zeros(size(div,1),size(div,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    fn = find(~isnan(Wdiv)); % find non NaN values
    div1sigma_20C(n,k) = std(Wdiv(fn)); % take the standard deviation of all non-NaN divergence values
    fd = find(Wdiv>0);	% find divergent points in space and time
    meanposdiv_20C(n,k) = sum((Wdiv(fd).*rep_cellarea(fd))/sum(rep_cellarea(fd))); % find the mean of the positive divergence, using the area-weighted array 
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    map_WI2O = squeeze(nansum(WI2O*30,3)); % 30 day months
    total_WI2O_20C(n,k) = squeeze(nansum(nansum(map_WI2O.*cellarea))); % take total I2O over all days in the winter
    ffrz = find(WI2O<0); % find freezing points in time and space
    WI2O_vol = WI2O.*cellarea*30; % 30 day months
    iceprod_20C(n,k) = squeeze(nansum(nansum(nansum(WI2O_vol(ffrz)))));
    M1 = zeros(size(WI2O)); M1(ffrz)=1;
    frz_area_days_20C(n,k) = squeeze(nansum(nansum(nansum(M1.*30.*cellarea))));
    fdf = intersect(ffrz,fd); % find divergent points that are also freezing points in time and space
    frz_pts_20C(n,k) = length(ffrz); % total number of freezing points in time and space through a winter
    iceprod_IF_20C(n,k) = squeeze(nansum(nansum(nansum(WI2O_vol(fdf))))); 
    Y = WI2O(fdf); X = [ones(length(Wdiv(fdf)),1),Wdiv(fdf)];
    [b,bint,r,rint,stats] = regress(Y,X);
    B0_20C(n,k) = b(1); B1_20C(n,k) = b(2); R2_20C(n,k) = stats(1);
    clearvars Wdiv WI2O
end
    fprintf('finished 20C ensemble # %d\n',k)
% now save the first part of the crossover year
n = 86;
    Wdiv = zeros(size(div,1),size(div,2),7);
    WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
     % take winter months
    Wdiv(:,:,1:3) = msk_div(:,:,12*(n-1)+10:12*(n-1)+12);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O(:,:,1:3) = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day

    clearvars msk_I2O msk_div I2O div 

% RCP8.5
if k >= 34
	div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-210012.nc'),'divu',[1,1,1],[291,48,900]);
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-210012.nc'),'fresh',[1,1,1],[291,48,900]);
else
	div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-208012.nc'),'divu');
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-208012.nc'),'fresh');
end
msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
for i = 1:size(div,3)
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end

% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
    Wdiv(:,:,4:7) = msk_div(:,:,1:4);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O(:,:,4:7) = msk_I2O(:,:,1:4).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    fn = find(~isnan(Wdiv)); % find non NaN values
    div1sigma_20C(n,k) = std(Wdiv(fn)); % take the standard deviation of all non-NaN divergence values
    fd = find(Wdiv>0);  % find divergent points in space and time
    meanposdiv_20C(n,k) = sum((Wdiv(fd).*rep_cellarea(fd))/sum(rep_cellarea(fd))); % find the mean of the positive divergence, using the area-weighted array 
    map_WI2O = squeeze(nansum(WI2O*30,3)); % 30 day months
    total_WI2O_20C(n,k) = squeeze(nansum(nansum(map_WI2O.*cellarea))); % take total I2O over all days in the winter
    ffrz = find(WI2O<0); % find freezing points in time and space
    WI2O_vol = WI2O.*cellarea*30; % 30 day months
    iceprod_20C(n,k) = squeeze(nansum(nansum(nansum(WI2O_vol(ffrz)))));
    M1 = zeros(size(WI2O)); M1(ffrz)=1;
    frz_area_days_20C(n,k) = squeeze(nansum(nansum(nansum(M1.*30.*cellarea))));
    fdf = intersect(ffrz,fd); % find divergent points that are also freezing points in time and space
    frz_pts_20C(n,k) = length(ffrz); % total number of freezing points in time and space through a winter
    iceprod_IF_20C(n,k) = squeeze(nansum(nansum(nansum(WI2O_vol(fdf)))));
    Y = WI2O(fdf); X = [ones(length(Wdiv(fdf)),1),Wdiv(fdf)];
    [b,bint,r,rint,stats] = regress(Y,X);
    B0_20C(n,k) = b(1); B1_20C(n,k) = b(2); R2_20C(n,k) = stats(1);
    clearvars Wdiv WI2O

% now for the rest of RCP8.5    
    
for n = 1:73 % no of years - 2 incomplete winter seasons
    Wdiv = zeros(size(div,1),size(div,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16); %.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    fn = find(~isnan(Wdiv)); % find non NaN values
    div1sigma_RCP85(n,k) = std(Wdiv(fn)); % take the standard deviation of all non-NaN divergence values
    fd = find(Wdiv>0);  % find divergent points in space and time
    meanposdiv_RCP85(n,k) = sum((Wdiv(fd).*rep_cellarea(fd))/sum(rep_cellarea(fd))); % find the mean of the positive divergence, using the area-weighted array 
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
    map_WI2O = squeeze(nansum(WI2O*30,3)); % 30 day months
    total_WI2O_RCP85(n,k) = squeeze(nansum(nansum(map_WI2O.*cellarea))); % take total I2O over all days in the winter
    ffrz = find(WI2O<0);  % find freezing points in time and space
    WI2O_vol = WI2O.*cellarea*30; % 30 day months
    iceprod_RCP85(n,k) = squeeze(nansum(nansum(nansum(WI2O_vol(ffrz)))));
    M1 = zeros(size(WI2O)); M1(ffrz)=1;
    frz_area_days_RCP85(n,k) = squeeze(nansum(nansum(nansum(M1.*30.*cellarea))));
    fdf = intersect(ffrz,fd); % find divergent points that are also freezing points in time and space
    frz_pts_RCP85(n,k) = length(ffrz); % total number of freezing points in time and space through a winter
    fdf = intersect(ffrz,fd); % find divergent points that are also freezing points in time and space
    iceprod_IF_RCP85(n,k) = squeeze(nansum(nansum(nansum(WI2O_vol(fdf))))); 
    Y = WI2O(fdf); X = [ones(length(Wdiv(fdf)),1),Wdiv(fdf)];  % construct regression
    [b,bint,r,rint,stats] = regress(Y,X);
    B0_RCP85(n,k) = b(1); B1_RCP85(n,k) = b(2); R2_RCP85(n,k) = stats(1);
    clearvars Wdiv WI2O
end
fprintf('finished RCP8.5 ensemble # %d\n',k)
end

save /home/ocean_personal_data/samc/CESM/variables/CESM_div_I2O_reg_20C_RCP85.mat meanposdiv_20C total_WI2O_20C frz_pts_20C B0_20C B1_20C R2_20C div1sigma_20C iceprod_20C iceprod_IF_20C frz_area_days_20C meanposdiv_RCP85 total_WI2O_RCP85 frz_pts_RCP85 B0_RCP85 B1_RCP85 R2_RCP85 div1sigma_RCP85 iceprod_RCP85 iceprod_IF_RCP85 frz_area_days_RCP85 
