% script to find changing winter cycle of sea ice freeze vs melt in KL seas--ensemble members 36 to 40
clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/hi/
addpath /home/ocean_shared_data2/CESM_arctic/aice/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_shared_data2/CESM_arctic/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/hs/

load KLmask.mat mask1xx mask
load fnames40.mat

dytxx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.fresh_nh.192001-200512.nc','dyt');
dxtxx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.fresh_nh.192001-200512.nc','dxt');
cellareaxx = dytxx.*dxtxx;
rep_cellareaxx = repmat(cellareaxx,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
rep_cellarea = repmat(cellarea,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise

monthly_I2O_20C = zeros(84,7,40);
monthly_melt_20C = zeros(84,7,40);
monthly_iceprod_20C = zeros(84,7,40);
monthly_I2O_RCP85 = zeros(84,7,40);
monthly_melt_RCP85 = zeros(84,7,40);
monthly_iceprod_RCP85 = zeros(84,7,40);
for k = 36:40    % total number of ensembles that we will use
% 20th CENTURY
       I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh');

       msk_I2O = zeros(size(I2O));
	for i = 1:size(I2O,3)
	msk_I2O(:,:,i) = I2O(:,:,i).*mask1xx;
	end

	for n = 1:85 % no of years - 1 incomplete winter seasons
        WI2O =zeros(size(I2O,1),size(I2O,2),7);
	% take winter months
        WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01*30;% assuming the I2O has ice conc built into it, converting cm /day --> m/day, and m/day --> m/month
	ffrz = find(WI2O<0); % find freezing points in time and space	
	fm = find(WI2O>0); % find freezing points in time and space	
	% have to replace the cells prior to processing as retrieving monthly output - must preserve array structure
	WI2O_frz = NaN*zeros(size(WI2O));
	WI2O_melt = NaN*zeros(size(WI2O));
	WI2O_frz(ffrz) = WI2O(ffrz);
	WI2O_melt(fm) = WI2O(fm);
        fn = find(~isnan(WI2O));

	monthly_I2O_20C(n,:,k) = squeeze(nansum(nansum(WI2O.*rep_cellareaxx)));
	monthly_iceprod_20C(n,:,k) = squeeze(nansum(nansum(WI2O_frz.*rep_cellareaxx)));
	monthly_melt_20C(n,:,k) = squeeze(nansum(nansum(WI2O_melt.*rep_cellareaxx)));
	end

 fprintf('finished 20C ensemble # %d\n',k)
n = 86;
WI2O = zeros(size(I2O,1),size(I2O,2),7);
% take winter months
WI2O(:,:,1:3) = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01*30;% assuming the I2O has ice conc built into it, converting cm /day --> m/day, and m/day --> m/month
 % RCP8.5
	I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-210012.nc'),'fresh',[1,1,1],[320,49,900]);

        msk_I2O = zeros(size(I2O));
	for i = 1:size(I2O,3)
	msk_I2O(:,:,i) = I2O(:,:,i).*mask1xx;
	end
% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
        WI2O(:,:,4:7) = msk_I2O(:,:,1:4).*0.01*30; % assuming the I2O has ice conc built into it, converting cm /day --> m/day
      ffrz = find(WI2O<0); % find freezing points in time and space
        fm = find(WI2O>0); % find freezing points in time and space
        % have to replace the cells prior to processing as retrieving monthly output - must preserve array structure
        WI2O_frz = NaN*zeros(size(WI2O));
        WI2O_melt = NaN*zeros(size(WI2O));
        WI2O_frz(ffrz) = WI2O(ffrz);
        WI2O_melt(fm) = WI2O(fm);
        fn = find(~isnan(WI2O));

        monthly_I2O_20C(n,:,k) = squeeze(nansum(nansum(WI2O.*rep_cellareaxx)));
        monthly_iceprod_20C(n,:,k) = squeeze(nansum(nansum(WI2O_frz.*rep_cellareaxx)));
        monthly_melt_20C(n,:,k) = squeeze(nansum(nansum(WI2O_melt.*rep_cellareaxx)));



	for n = 1:73 % no of years - 2 incomplete winter seasons
        WI2O =zeros(size(I2O,1),size(I2O,2),7);
	% take winter months
        WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01*30;% assuming the I2O has ice conc built into it, converting cm /day --> m/day, and m/day --> m/month
	ffrz = find(WI2O<0); % find freezing points in time and space	
	fm = find(WI2O>0); % find freezing points in time and space	
	% have to replace the cells prior to processing as retrieving monthly output - must preverse array structure
	WI2O_frz = NaN*zeros(size(WI2O));
	WI2O_melt = NaN*zeros(size(WI2O));
	WI2O_frz(ffrz) = WI2O(ffrz);
	WI2O_melt(fm) = WI2O(fm);
        fn = find(~isnan(WI2O));

	monthly_I2O_RCP85(n,:,k) = squeeze(nansum(nansum(WI2O.*rep_cellareaxx)));
	monthly_iceprod_RCP85(n,:,k) = squeeze(nansum(nansum(WI2O_frz.*rep_cellareaxx)));
	monthly_melt_RCP85(n,:,k) = squeeze(nansum(nansum(WI2O_melt.*rep_cellareaxx)));
	end

 fprintf('finished RCP85 ensemble # %d\n',k)
 end

save /home/ocean_personal_data/samc/CESM/variables/CESM_monthly_I2O_36to40.mat monthly_I2O_20C monthly_iceprod_20C monthly_melt_20C...
monthly_I2O_RCP85 monthly_iceprod_RCP85 monthly_melt_RCP85
