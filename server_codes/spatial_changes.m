% spatial changes - taking a mean through the ensembles and through winter months.

clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/hi/
addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/aice/
addpath /home/ocean_shared_data2/CESM_arctic/hs/
addpath /home/ocean_shared_data2/CESM_arctic/Tair/
addpath /home/ocean_shared_data2/CESM_arctic/TEMP/
load KLmask_pop.mat; mask_pop = mask; clearvars mask;
load KLmask.mat
load fnames40.mat
load ArcticMask.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
rep_cellarea = repmat(cellarea,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise


for k = 1:35    % total number of ensembles that we will use
% 20th CENTURY
if k == 1
        div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.185001-200512.nc'),'divu',[1,1,841],[291,48,1032]);
        I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.185001-200512.nc'),'fresh',[1,1,841],[291,48,1032]);
        SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.185001-200512.nc'),'aice',[1,1,841],[291,48,1032])/100;
	SIT = ncread(strcat(fnames_20C{k},'.cice.h.hi_nh.185001-200512.nc'),'hi',[1,1,841],[291,48,1032]);
        snow = ncread(strcat(fnames_20C{k},'.cice.h.hs_nh.185001-200512.nc'),'hs',[1,1,841],[291,48,1032]);
        SAT = ncread(strcat(fnames_20C{k},'.cice.h.Tair_nh.185001-200512.nc'),'Tair',[1,1,841],[291,48,1032]);
	T10m = squeeze(ncread(strcat(fnames_20C{k},'.pop.h.TEMP.185001-200512.nc'),'TEMP',[1,1,1,841],[320,49,1,1032]));
else 	
        div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
        I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh');
        SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.192001-200512.nc'),'aice')/100;
        SIT = ncread(strcat(fnames_20C{k},'.cice.h.hi_nh.192001-200512.nc'),'hi');
        snow = ncread(strcat(fnames_20C{k},'.cice.h.hs_nh.192001-200512.nc'),'hs');
        SAT = ncread(strcat(fnames_20C{k},'.cice.h.Tair_nh.192001-200512.nc'),'Tair');
	T10m = squeeze(ncread(strcat(fnames_20C{k},'.pop.h.TEMP.192001-200512.nc'),'TEMP',[1,1,1,1],[320,49,1,1032]));
end
	 msk_I2O = zeros(size(I2O));
	 Arcticmsk_I2O = zeros(size(I2O));
	 for i = 1:size(div,3)
	 msk_I2O(:,:,i) = I2O(:,:,i).*mask;
	 Arcticmsk_I2O(:,:,i) = I2O(:,:,i).*ArcticMask; 
	 end

	 for n = 1:85 % no. of years minus 1 incomplete winter season
                 Wdiv =zeros(size(div,1),size(div,2),7);
                 WI2O =zeros(size(I2O,1),size(I2O,2),7);
                 WSIC =zeros(size(SIC,1),size(SIC,2),7);
                 WSIT =zeros(size(SIT,1),size(SIT,2),7);
                 Wsnow = zeros(size(snow,1),size(snow,2),7);
                 WSAT =zeros(size(SAT,1),size(SAT,2),7);
                 WSIA =zeros(size(SIC,1),size(SIC,2),7);
		 WT10m = zeros(size(T10m,1),size(T10m,2),7);	
	if k == 1 && n == 1	
	 MKL = zeros(size(SIC,1),size(SIC,2));
	 fnKL = find(~isnan(msk_I2O(:,:,1)));
	 MKL(fnKL) = cellarea(fnKL);
	 areaKL = sum(MKL(:));
	 MArctic = zeros(size(SIC,1),size(SIC,2));
	 fnArctic = find(~isnan(Arcticmsk_I2O(:,:,1)));
	 MArctic(fnArctic) = cellarea(fnArctic);
	 areaArctic = sum(MArctic(:));
	end

         % take winter months
         Wdiv = div(:,:,12*(n-1)+10:12*(n-1)+16);
	 Wdiv(Wdiv==0) = NaN; % Nan-ing non ice areas
         WI2O = I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_Arctic = Arcticmsk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_KL = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WSIC = SIC(:,:,12*(n-1)+10:12*(n-1)+16);
         SepSIC = SIC(:,:,12*(n-1)+9);
         WSIT = SIT(:,:,12*(n-1)+10:12*(n-1)+16);
         SepSIT = SIT(:,:,12*(n-1)+9);
         Wsnow = snow(:,:,12*(n-1)+10:12*(n-1)+16);
	 Wsnow(Wsnow==0) = NaN; % Nan-ing non ice areas
         WSAT = SAT(:,:,12*(n-1)+10:12*(n-1)+16);
	 SepT10m = T10m(:,:,12*(n-1)+9);
	 ffrz_A = find(WI2O_Arctic<0); % find freezing points in time and space
	 ffrz_KL = find(WI2O_KL<0); % find freezing points in time and space
	 iceprod_Arctic_20C(n,k) = sum(WI2O_Arctic(ffrz_A).*rep_cellarea(ffrz_A)*30);
	 iceprod_KL_20C(n,k) = sum(WI2O_KL(ffrz_KL).*rep_cellarea(ffrz_KL)*30);
	 ffrz = find(WI2O<0);
	 WI2O_frz = NaN*zeros(size(WI2O));
	 WI2O_frz(ffrz) = WI2O(ffrz);
	 
	 % take winter mean. Be aware of NaNs...
	 div_map_20C(:,:,n,k) = nanmean(Wdiv,3);
	 I2O_map_20C(:,:,n,k) = mean(WI2O,3);
	 I2O_frz_map_20C(:,:,n,k) = nanmean(WI2O_frz,3);
	 SIC_map_20C(:,:,n,k) = mean(WSIC,3);
	 SIT_map_20C(:,:,n,k) = mean(WSIT,3);
	 snow_map_20C(:,:,n,k) = nanmean(Wsnow,3);
	 SepSIC_map_20C(:,:,n,k) = SepSIC;
	 SepSIT_map_20C(:,:,n,k) = SepSIT;
	 SAT_map_20C(:,:,n,k) = mean(WSAT,3);
	 SepT10m_map_20C(:,:,n,k) = SepT10m;
 end
 fprintf('finished 20C ensemble # %d\n',k)

% now save the first part of the crossover year
n = 86;
 		 Wdiv =zeros(size(div,1),size(div,2),7);
                 WI2O =zeros(size(I2O,1),size(I2O,2),7);
                 WSIC =zeros(size(SIC,1),size(SIC,2),7);
                 WSIT =zeros(size(SIT,1),size(SIT,2),7);
                 Wsnow = zeros(size(snow,1),size(snow,2),7);
                 WSAT =zeros(size(SAT,1),size(SAT,2),7);
                 WSIA =zeros(size(SIC,1),size(SIC,2),7);
                 WT10m = zeros(size(T10m,1),size(T10m,2),7);
		 % take winter months
 	 Wdiv(:,:,1:3) = div(:,:,12*(n-1)+10:12*(n-1)+12);
         WI2O(:,:,1:3) = I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_Arctic(:,:,1:3) = Arcticmsk_I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_KL(:,:,1:3) = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WSIC(:,:,1:3) = SIC(:,:,12*(n-1)+10:12*(n-1)+12);
         SepSIC = SIC(:,:,12*(n-1)+9);
         WSIT(:,:,1:3) = SIT(:,:,12*(n-1)+10:12*(n-1)+12);
         SepSIT = SIT(:,:,12*(n-1)+9);
         Wsnow(:,:,1:3) = snow(:,:,12*(n-1)+10:12*(n-1)+12);
         WSAT(:,:,1:3) = SAT(:,:,12*(n-1)+10:12*(n-1)+12);
         SepT10m = T10m(:,:,12*(n-1)+9);
 % RCP8.5
if k >= 34
        div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-210012.nc'),'divu',[1,1,1],[291,48,900]);
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-210012.nc'),'fresh',[1,1,1],[291,48,900]);
        SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-210012.nc'),'aice',[1,1,1],[291,48,900])/100;
	SIT = ncread(strcat(fnames_RCP85{k},'.cice.h.hi_nh.200601-210012.nc'),'hi',[1,1,1],[291,48,900]);
        snow = ncread(strcat(fnames_RCP85{k},'.cice.h.hs_nh.200601-210012.nc'),'hs',[1,1,1],[291,48,900]);
        SAT = ncread(strcat(fnames_RCP85{k},'.cice.h.Tair_nh.200601-210012.nc'),'Tair',[1,1,1],[291,48,900]);
        T10m = squeeze(ncread(strcat(fnames_RCP85{k},'.pop.h.TEMP.200601-210012.nc'),'TEMP',[1,1,1,1],[320,49,1,900]));

else
         div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-208012.nc'),'divu');
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-208012.nc'),'fresh');
        SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-208012.nc'),'aice')/100;
        SIT = ncread(strcat(fnames_RCP85{k},'.cice.h.hi_nh.200601-208012.nc'),'hi');
        snow = ncread(strcat(fnames_RCP85{k},'.cice.h.hs_nh.200601-208012.nc'),'hs');
        SAT = ncread(strcat(fnames_RCP85{k},'.cice.h.Tair_nh.200601-208012.nc'),'Tair');
	T10m = squeeze(ncread(strcat(fnames_RCP85{k},'.pop.h.TEMP.200601-208012.nc'),'TEMP',[1,1,1,1],[320,49,1,900]));
end

         msk_I2O = zeros(size(I2O));
         Arcticmsk_I2O = zeros(size(I2O));
         for i = 1:size(div,3)
         msk_I2O(:,:,i) = I2O(:,:,i).*mask;
         Arcticmsk_I2O(:,:,i) = I2O(:,:,i).*ArcticMask;
         end   

% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
	 Wdiv(:,:,4:7) = div(:,:,1:4);
         WI2O(:,:,4:7) = I2O(:,:,1:4).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_Arctic(:,:,4:7) = Arcticmsk_I2O(:,:,1:4).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_KL(:,:,4:7) = msk_I2O(:,:,1:4).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WSIC(:,:,4:7) = SIC(:,:,1:4);
         WSIT(:,:,4:7) = SIT(:,:,1:4);
         Wsnow(:,:,4:7) = snow(:,:,1:4);
         WSAT(:,:,4:7) = SAT(:,:,1:4);
	 Wdiv(Wdiv==0) = NaN; % Nan-ing non ice areas
	 Wsnow(Wsnow==0) = NaN; % Nan-ing non ice areas
	 ffrz_A = find(WI2O_Arctic<0); % find freezing points in time and space
         ffrz_KL = find(WI2O_KL<0); % find freezing points in time and space
         iceprod_Arctic_20C(n,k) = sum(WI2O_Arctic(ffrz_A).*rep_cellarea(ffrz_A)*30);
         iceprod_KL_20C(n,k) = sum(WI2O_KL(ffrz_KL).*rep_cellarea(ffrz_KL)*30);
         ffrz = find(WI2O<0);
         WI2O_frz = NaN*zeros(size(WI2O));
         WI2O_frz(ffrz) = WI2O(ffrz);

         % take winter mean. Be aware of NaNs...
         div_map_20C(:,:,n,k) = nanmean(Wdiv,3);
         I2O_map_20C(:,:,n,k) = mean(WI2O,3);
         I2O_frz_map_20C(:,:,n,k) = nanmean(WI2O_frz,3);
         SIC_map_20C(:,:,n,k) = mean(WSIC,3);
         SIT_map_20C(:,:,n,k) = mean(WSIT,3);
         snow_map_20C(:,:,n,k) = nanmean(Wsnow,3);
         SepSIC_map_20C(:,:,n,k) = SepSIC;
         SepSIT_map_20C(:,:,n,k) = SepSIT;
         SAT_map_20C(:,:,n,k) = mean(WSAT,3);
         SepT10m_map_20C(:,:,n,k) = SepT10m;

	for n = 1:73 % no. of years mins 2 incomplete winter seasons
       	         Wdiv =zeros(size(div,1),size(div,2),7);
                 WI2O =zeros(size(I2O,1),size(I2O,2),7);
                 WSIC =zeros(size(SIC,1),size(SIC,2),7);
                 WSIT =zeros(size(SIT,1),size(SIT,2),7);
                 Wsnow = zeros(size(snow,1),size(snow,2),7);
                 WSAT =zeros(size(SAT,1),size(SAT,2),7);
                 WSIA =zeros(size(SIC,1),size(SIC,2),7);
		 WT10m = zeros(size(T10m,1),size(T10m,2),7);	

         % take winter months
         Wdiv = div(:,:,12*(n-1)+10:12*(n-1)+16);
	 Wdiv(Wdiv==0) = NaN; % Nan-ing non ice areas
         WI2O = I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_Arctic = Arcticmsk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WI2O_KL = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
         WSIC = SIC(:,:,12*(n-1)+10:12*(n-1)+16);
         SepSIC = SIC(:,:,12*(n-1)+9);
         WSIT = SIT(:,:,12*(n-1)+10:12*(n-1)+16);
         SepSIT = SIT(:,:,12*(n-1)+9);
         Wsnow = snow(:,:,12*(n-1)+10:12*(n-1)+16);
	 Wsnow(Wsnow==0) = NaN; % Nan-ing non ice areas
         WSAT = SAT(:,:,12*(n-1)+10:12*(n-1)+16);
	 SepT10m = T10m(:,:,12*(n-1)+9);
	 ffrz_A = find(WI2O_Arctic<0); % find freezing points in time and space
         ffrz_KL = find(WI2O_KL<0); % find freezing points in time and space
         iceprod_Arctic_RCP85(n,k) = sum(WI2O_Arctic(ffrz_A).*rep_cellarea(ffrz_A)*30);
         iceprod_KL_RCP85(n,k) = sum(WI2O_KL(ffrz_KL).*rep_cellarea(ffrz_KL)*30);
	 ffrz = find(WI2O<0);
	 WI2O_frz = NaN*zeros(size(WI2O));
	 WI2O_frz(ffrz) = WI2O(ffrz);

         % take winter mean. Be aware of NaNs...
         div_map_RCP85(:,:,n,k) = nanmean(Wdiv,3);
         I2O_map_RCP85(:,:,n,k) = mean(WI2O,3);
	 I2O_frz_map_RCP85(:,:,n,k) = nanmean(WI2O_frz,3);
         SIC_map_RCP85(:,:,n,k) = mean(WSIC,3);
         SIT_map_RCP85(:,:,n,k) = mean(WSIT,3);
         snow_map_RCP85(:,:,n,k) = nanmean(Wsnow,3);
         SepSIC_map_RCP85(:,:,n,k) = SepSIC;
         SepSIT_map_RCP85(:,:,n,k) = SepSIT;
	 SAT_map_RCP85(:,:,n,k) = mean(WSAT,3);
	 SepT10m_map_RCP85(:,:,n,k) = SepT10m;
 end
 fprintf('finished RCP85 ensemble # %d\n',k)
end

% ensemble means

	div_map_20C = mean(div_map_20C,4);
	I2O_map_20C = mean(I2O_map_20C,4);
	SIC_map_20C = mean(SIC_map_20C,4);
	SIT_map_20C = mean(SIT_map_20C,4);
	snow_map_20C = mean(snow_map_20C,4);
	SepSIC_map_20C = mean(SepSIC_map_20C,4);
	SepSIT_map_20C = mean(SepSIT_map_20C,4);
	SAT_map_20C = mean(SAT_map_20C,4);
	I2O_frz_map_20C = nanmean(I2O_frz_map_20C,4);
	SepT10m_map_20C = mean(SepT10m_map_20C,4);

	div_map_RCP85 = mean(div_map_RCP85,4);
	I2O_map_RCP85 = mean(I2O_map_RCP85,4);
	SIC_map_RCP85 = mean(SIC_map_RCP85,4);
	SIT_map_RCP85 = mean(SIT_map_RCP85,4);
	snow_map_RCP85 = mean(snow_map_RCP85,4);
	SepSIC_map_RCP85 = mean(SepSIC_map_RCP85,4);
	SepSIT_map_RCP85 = mean(SepSIT_map_RCP85,4);
	SAT_map_RCP85 = mean(SAT_map_RCP85,4);
	I2O_frz_map_RCP85 = nanmean(I2O_frz_map_RCP85,4);
	SepT10m_map_RCP85 = mean(SepT10m_map_RCP85,4);

save /home/ocean_personal_data/samc/CESM/variables/CESM_spatial_maps.mat div_map_20C I2O_map_20C SIC_map_20C SIT_map_20C snow_map_20C SepSIC_map_20C SepSIT_map_20C div_map_RCP85 I2O_map_RCP85 SIC_map_RCP85 SIT_map_RCP85 snow_map_RCP85 SepSIC_map_RCP85 SepSIT_map_RCP85 SAT_map_20C SAT_map_RCP85 I2O_frz_map_20C I2O_frz_map_RCP85 iceprod_Arctic_20C iceprod_Arctic_RCP85 iceprod_KL_20C iceprod_KL_RCP85 areaKL areaArctic SepT10m_map_20C SepT10m_map_RCP85
