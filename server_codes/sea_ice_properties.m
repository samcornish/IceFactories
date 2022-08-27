% script to find changing winter cycle of sea ice area in KL seas
clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/hi/
addpath /home/ocean_shared_data2/CESM_arctic/aice/
addpath /home/ocean_shared_data2/CESM_arctic/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/hs/

load KLmask.mat
load fnames40.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
monthly_SIT_20C = zeros(84,8,35);
monthly_SIA_20C = zeros(84,8,35); % initialise
monthly_SIT_RCP85 = zeros(73,8,35);
monthly_SIA_RCP85 = zeros(73,8,35); % initialise
for k = 1:35	% total number of ensembles that we will use
% 20th CENTURY
if k == 1
	SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.185001-200512.nc'),'aice',[1,1,841],[291,48,1032])/100;
	SIT = ncread(strcat(fnames_20C{k},'.cice.h.hi_nh.185001-200512.nc'),'hi',[1,1,841],[291,48,1032]);
	snow = ncread(strcat(fnames_20C{k},'.cice.h.hs_nh.185001-200512.nc'),'hs',[1,1,841],[291,48,1032]);
else
	SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.192001-200512.nc'),'aice')/100;
	SIT = ncread(strcat(fnames_20C{k},'.cice.h.hi_nh.192001-200512.nc'),'hi');
	snow = ncread(strcat(fnames_20C{k},'.cice.h.hs_nh.192001-200512.nc'),'hs');
end
msk_SIT = zeros(size(SIT));	
msk_SIC = zeros(size(SIC));
msk_snow = zeros(size(snow));
for i = 1:size(SIC,3)
msk_SIT(:,:,i) = SIT(:,:,i).*mask;
msk_SIC(:,:,i) = SIC(:,:,i).*mask;
msk_snow(:,:,i) = snow(:,:,i).*mask;
end

for n = 1:85 % no of years - 2 incomplete winter seasons
    WSIT = NaN*zeros(size(SIT,1),size(SIT,2),8);
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),8);
    WSIA = NaN*zeros(size(SIC,1),size(SIC,2),8);
    Wsnow = NaN*zeros(size(snow,1),size(snow,2),8);
    % take winter months
    WSIT = msk_SIT(:,:,12*(n-1)+9:12*(n-1)+16);
    WSIC = msk_SIC(:,:,12*(n-1)+9:12*(n-1)+16);
    Wsnow = msk_snow(:,:,12*(n-1)+9:12*(n-1)+16);
    rep_cellarea = repmat(cellarea,[1,1,8]);
    WSIA = WSIC.*rep_cellarea;
    fp = find(WSIC>0.15); % find cells with more than 15% sea ice
%	WSIT_wt_conc = WSIT.*WSIC;  % not clear this is necessary... It is not: 'hi' gives the ice volume/unit grid cell area
	M = ones(size(WSIT));
	fn = find(isnan(WSIT));
	M(fn) = NaN; 
	area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));	% sum of all areas where SIT is recorded
    monthly_SIT_20C(n,:,k) = squeeze(nansum(nansum(WSIT.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    monthly_SIV_20C(n,:,k) = squeeze(nansum(nansum(WSIT.*rep_cellarea))); 
    monthly_SV_20C(n,:,k) =  squeeze(nansum(nansum(Wsnow.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    monthly_SIA_20C(n,:,k) = squeeze(nansum(nansum(WSIA)));
    monthly_hs_20C(n,:,k) = squeeze(nansum(nansum(((Wsnow.*rep_cellarea)./WSIA).*rep_cellarea)))./area_sum;	% should be the mean snow depth on ice, not volume/cell area
    Wh_onice = (Wsnow.*rep_cellarea)./(WSIC.*rep_cellarea);
    snow_over_ice_extent_20C(n,k) = sum((Wh_onice(fp).*rep_cellarea(fp))/sum(rep_cellarea(fp))); 
end   

fprintf('finished 20C ensemble # %d\n',k)
% now save the first part of the crossover year
n = 86;
    WSIT = NaN*zeros(size(SIT,1),size(SIT,2),8);
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),8);
    WSIA = NaN*zeros(size(SIC,1),size(SIC,2),8);
    Wsnow = NaN*zeros(size(snow,1),size(snow,2),8);
% take winter months
    WSIT(:,:,1:4) = msk_SIT(:,:,12*(n-1)+9:12*(n-1)+12);
    WSIC(:,:,1:4) = msk_SIC(:,:,12*(n-1)+9:12*(n-1)+12);
    Wsnow(:,:,1:4) = msk_snow(:,:,12*(n-1)+9:12*(n-1)+12);


% RCP8.5
if k >= 34
	SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-210012.nc'),'aice',[1,1,1],[291,48,900])/100;
        SIT = ncread(strcat(fnames_RCP85{k},'.cice.h.hi_nh.200601-210012.nc'),'hi',[1,1,1],[291,48,900]);
        snow = ncread(strcat(fnames_RCP85{k},'.cice.h.hs_nh.200601-210012.nc'),'hs',[1,1,1],[291,48,900]);
else
	SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-208012.nc'),'aice')/100;
	SIT = ncread(strcat(fnames_RCP85{k},'.cice.h.hi_nh.200601-208012.nc'),'hi');
	snow = ncread(strcat(fnames_RCP85{k},'.cice.h.hs_nh.200601-208012.nc'),'hs');
end
msk_SIT = zeros(size(SIT));	
msk_SIC = zeros(size(SIC));
msk_snow = zeros(size(snow));
for i = 1:size(SIC,3)
msk_SIT(:,:,i) = SIT(:,:,i).*mask;
msk_SIC(:,:,i) = SIC(:,:,i).*mask;
msk_snow(:,:,i) = snow(:,:,i).*mask;
end

% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
    WSIT(:,:,5:8) = msk_SIT(:,:,1:4);
    WSIC(:,:,5:8) = msk_SIC(:,:,1:4);
    Wsnow(:,:,5:8) = msk_snow(:,:,1:4);
    WSIA = WSIC.*rep_cellarea;
    fp = find(WSIC>0.15); % find cells with more than 15% sea ice
%       WSIT_wt_conc = WSIT.*WSIC;  % not clear this is necessary... It is not: 'hi' gives the ice volume/unit grid cell area
        M = ones(size(WSIT));
        fn = find(isnan(WSIT));
        M(fn) = NaN;
        area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));    % sum of all areas where SIT is recorded
    monthly_SIT_20C(n,:,k) = squeeze(nansum(nansum(WSIT.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    monthly_SIV_20C(n,:,k) = squeeze(nansum(nansum(WSIT.*rep_cellarea)));
    monthly_SV_20C(n,:,k) =  squeeze(nansum(nansum(Wsnow.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    monthly_SIA_20C(n,:,k) = squeeze(nansum(nansum(WSIA)));
    monthly_hs_20C(n,:,k) = squeeze(nansum(nansum(((Wsnow.*rep_cellarea)./WSIA).*rep_cellarea)))./area_sum;     % should be the mean snow depth on ice, not volume/cell area
    Wh_onice = (Wsnow.*rep_cellarea)./(WSIC.*rep_cellarea);
    snow_over_ice_extent_20C(n,k) = sum((Wh_onice(fp).*rep_cellarea(fp))/sum(rep_cellarea(fp)));


for n = 1:73 % no of years - 2 incomplete winter seasons
    WSIT = NaN*zeros(size(SIT,1),size(SIT,2),8);
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),8);
    WSIA = NaN*zeros(size(SIC,1),size(SIC,2),8);
    Wsnow = NaN*zeros(size(snow,1),size(snow,2),8);
    % take winter months
    WSIT = msk_SIT(:,:,12*(n-1)+9:12*(n-1)+16);
    WSIC = msk_SIC(:,:,12*(n-1)+9:12*(n-1)+16);
    Wsnow = msk_snow(:,:,12*(n-1)+9:12*(n-1)+16);
    WSIA = WSIC.*rep_cellarea;
    fp = find(WSIC>0.15); % find cells with more than 15% sea ice
%	WSIT_wt_conc = WSIT.*WSIC;  % not clear this is necessary... It is not: 'hi' gives the ice volume/unit grid cell area
	M = ones(size(WSIT));
	fn = find(isnan(WSIT));
	M(fn) = NaN; 
	area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));	% sum of all areas where SIT is recorded
    monthly_SIT_RCP85(n,:,k) = squeeze(nansum(nansum(WSIT.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    monthly_SIV_RCP85(n,:,k) = squeeze(nansum(nansum(WSIT.*rep_cellarea))); 
    monthly_SV_RCP85(n,:,k) =  squeeze(nansum(nansum(Wsnow.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    monthly_SIA_RCP85(n,:,k) = squeeze(nansum(nansum(WSIA)));
    monthly_hs_RCP85(n,:,k) = squeeze(nansum(nansum(((Wsnow.*rep_cellarea)./WSIA).*rep_cellarea)))./area_sum;	% should be the mean snow depth on ice, not volume/cell area
    Wh_onice = (Wsnow.*rep_cellarea)./(WSIC.*rep_cellarea);
    snow_over_ice_extent_RCP85(n,k) = sum((Wh_onice(fp).*rep_cellarea(fp))/sum(rep_cellarea(fp))); 
end

fprintf('finished RCP8.5 ensemble # %d\n',k)
end
total_area = area_sum(1);
save /home/ocean_personal_data/samc/CESM/variables/CESM_SIT_SIC_SIV_20C_RCP85_Sep.mat monthly_SIA_20C monthly_SIT_20C monthly_SIV_20C monthly_SV_20C monthly_SIA_RCP85 monthly_SIT_RCP85 monthly_SIV_RCP85 monthly_SV_RCP85 monthly_hs_20C monthly_hs_RCP85 total_area snow_over_ice_extent_20C snow_over_ice_extent_RCP85
