% master script to collect all params relating to ice growth rate

clear 
close all

addpath /home/ocean_shared_data2/CESM_arctic/hi/
addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/aice/
addpath /home/ocean_shared_data2/CESM_arctic/hs/
addpath /home/ocean_shared_data2/CESM_arctic/Tair/
addpath /home/ocean_shared_data2/CESM_arctic/fhocn_ai/

load KLmask.mat mask1xx mask
load fnames40.mat

dytxx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.divu_nh.192001-200512.nc','dyt');
dxtxx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.divu_nh.192001-200512.nc','dxt');
cellareaxx = dytxx.*dxtxx;
rep_cellareaxx = repmat(cellareaxx,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise
dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.hi_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.hi_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
rep_cellarea = repmat(cellarea,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise

for k = 36:40    % total number of ensembles that we will use
% 20th CENTURY
        div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
        I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh');
        SIC = ncread(strcat(fnames_20C{k},'.cice.h.aice_nh.192001-200512.nc'),'aice')/100;
        snow = ncread(strcat(fnames_20C{k},'.cice.h.hs_nh.192001-200512.nc'),'hs');
	SAT = ncread(strcat(fnames_20C{k},'.cice.h.Tair_nh.192001-200512.nc'),'Tair');
	Fw = ncread(strcat(fnames_20C{k},'.cice.h.fhocn_ai_nh.192001-200512.nc'),'fhocn_ai');%W/m^2

	msk_div = zeros(size(div));
	msk_I2O = zeros(size(I2O));
	msk_SIC = zeros(size(SIC));
	msk_snow = zeros(size(snow));
	msk_SAT = zeros(size(SAT));
	msk_Fw = zeros(size(Fw));
	for i = 1:size(div,3)
	msk_div(:,:,i) = div(:,:,i).*mask1xx;
	msk_I2O(:,:,i) = I2O(:,:,i).*mask1xx;
	msk_SIC(:,:,i) = SIC(:,:,i).*mask1xx;
	msk_snow(:,:,i) = snow(:,:,i).*mask1xx;
	msk_SAT(:,:,i) = SAT(:,:,i).*mask1xx;
	msk_Fw(:,:,i) = Fw(:,:,i).*mask1xx;
	end

	for n = 1:85 % no. of years minus 1 incomplete winter season
		 Wdiv = NaN*zeros(size(div,1),size(div,2),7);
		 WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
		 WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
		 SepSIC = NaN*zeros(size(SIC,1),size(SIC,2)); 
		 Wsnow = NaN*zeros(size(snow,1),size(snow,2),7);
		 WSAT = NaN*zeros(size(SAT,1),size(SAT,2),7);
		 WFw = NaN*zeros(size(SIC,1),size(SIC,2),7);
	 % take winter months
	 Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);
	 WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
	 WSIC = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+16);
	 SepSIC = msk_SIC(:,:,12*(n-1)+9);
	 Wsnow = msk_snow(:,:,12*(n-1)+10:12*(n-1)+16);
	 WSAT = msk_SAT(:,:,12*(n-1)+10:12*(n-1)+16);
	 WFw = msk_Fw(:,:,12*(n-1)+10:12*(n-1)+16);
	 
	 fn = find(~isnan(WI2O)); % find non NaN values; this finds values over the ocean - however it looks like SAT is also NaN over land as from cice
	 fp = find(WSIC>0.15); % find cells with more than 15% sea ice
	 fnd = find(~isnan(Wdiv)); % find non NaN values of divergence
	 ffrz = find(WI2O<0); % find freezing points in time and space
	 fd = find(Wdiv>0);  % find divergent points in space and time
	 fpd = intersect(fd,fp); % find cells where >15% ice and divergent
	 fdf = intersect(fd,ffrz); % freezing and divergent points
	 fpf = intersect(fp,ffrz); % points where ice >15% and freezing
        % area-weighted growth rate only over cells where ice is growing
         GR_20C(n,k) = sum((WI2O(ffrz).*rep_cellareaxx(ffrz))/sum(rep_cellareaxx(ffrz)));
        % total ice production (m^3). Need factor 30 to convert from mean daily to mean monthly. 
	 iceprod_20C(n,k) = sum(WI2O(ffrz).*rep_cellareaxx(ffrz))*30;
	 % total positive area diverged each winter in m^2
	 pos_areadiv_20C(n,k) = sum((Wdiv(fd)*(30/100)).*rep_cellareaxx(fd)); % factor (30/100) converts %/day to /month
	 % net area diverged each winter in m^2
	 net_areadiv_20C(n,k) = sum((Wdiv(fnd)*(30/100)).*rep_cellareaxx(fnd)); % factor (30/100) converts %/day to /month
         % snow - count over >15% ice cells
         hs_20C(n,k) = sum((Wsnow(fp).*rep_cellareaxx(fp))/sum(rep_cellareaxx(fp)));
	 % September sea ice area
	 SepSIA_20C(n,k) = squeeze(nansum(nansum(SepSIC.*cellareaxx)));
         % temperature - count over all ocean cells
	 SAT_20C(n,k) = sum((WSAT(fn).*rep_cellareaxx(fn))/sum(rep_cellareaxx(fn)));
	 % freezing area days
	 M1 = ones(size(WI2O))*30;
	 frz_area_days_20C(n,k) = sum(M1(ffrz).*rep_cellareaxx(ffrz));
	 % heat fluxes - count and normalise over days where ice present at >15%
	 Fw_20C(n,k) = sum((WFw(ffrz).*rep_cellareaxx(ffrz))/sum(rep_cellareaxx(ffrz)));
	 end
 fprintf('finished 20C ensemble # %d\n',k)
% now save the first part of the crossover year
n = 86;
		Wdiv = NaN*zeros(size(div,1),size(div,2),7);
                WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
		WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
		SepSIC = NaN*zeros(size(SIC,1),size(SIC,2)); 
                Wsnow = NaN*zeros(size(snow,1),size(snow,2),7);
                WSAT = NaN*zeros(size(SAT,1),size(SAT,2),7);
                WFw = NaN*zeros(size(SIC,1),size(SIC,2),7);
% take winter months
         Wdiv(:,:,1:3) = msk_div(:,:,12*(n-1)+10:12*(n-1)+12);
         WI2O(:,:,1:3) = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+12).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
	 WSIC(:,:,1:3) = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+12); 
         SepSIC(:,:) = msk_SIC(:,:,12*(n-1)+9);
	 Wsnow(:,:,1:3) = msk_snow(:,:,12*(n-1)+10:12*(n-1)+12);
         WSAT(:,:,1:3) = msk_SAT(:,:,12*(n-1)+10:12*(n-1)+12);
         WFw(:,:,1:3) = msk_Fw(:,:,12*(n-1)+10:12*(n-1)+12);

% RCP8.5
        div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-210012.nc'),'divu',[1,1,1],[320,49,900]);
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-210012.nc'),'fresh',[1,1,1],[320,49,900]);
        SIC = ncread(strcat(fnames_RCP85{k},'.cice.h.aice_nh.200601-210012.nc'),'aice',[1,1,1],[320,49,900])/100;
        snow = ncread(strcat(fnames_RCP85{k},'.cice.h.hs_nh.200601-210012.nc'),'hs',[1,1,1],[320,49,900]);
        SAT = ncread(strcat(fnames_RCP85{k},'.cice.h.Tair_nh.200601-210012.nc'),'Tair',[1,1,1],[320,49,900]);
        Fw = ncread(strcat(fnames_RCP85{k},'.cice.h.fhocn_ai_nh.200601-210012.nc'),'fhocn_ai',[1,1,1],[320,49,900]);%W/m^2
	
        msk_div = zeros(size(div));
        msk_I2O = zeros(size(I2O));
        msk_SIC = zeros(size(SIC));
        msk_snow = zeros(size(snow));
        msk_SAT = zeros(size(SAT));
        msk_Fw = zeros(size(Fw));
        for i = 1:size(div,3)
        msk_div(:,:,i) = div(:,:,i).*mask1xx;
        msk_I2O(:,:,i) = I2O(:,:,i).*mask1xx;
        msk_SIC(:,:,i) = SIC(:,:,i).*mask1xx;
        msk_snow(:,:,i) = snow(:,:,i).*mask1xx;
        msk_SAT(:,:,i) = SAT(:,:,i).*mask1xx;
        msk_Fw(:,:,i) = Fw(:,:,i).*mask1xx;
        end
% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
	 Wdiv(:,:,4:7) = msk_div(:,:,1:4);
         WI2O(:,:,4:7) = msk_I2O(:,:,1:4).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
	 WSIC(:,:,4:7) = msk_SIC(:,:,1:4);
         Wsnow(:,:,4:7) = msk_snow(:,:,1:4);
         WSAT(:,:,4:7) = msk_SAT(:,:,1:4);
         WFw(:,:,4:7) = msk_Fw(:,:,1:4);

	 fn = find(~isnan(WI2O)); % find non NaN values; this finds values over the ocean - however it looks like SAT is also NaN over land as from cice
	 fp = find(WSIC>0.15); % find cells with more than 15% sea ice
	 fnd = find(~isnan(Wdiv)); % find non NaN values of divergence
         ffrz = find(WI2O<0); % find freezing points in time and space
         fd = find(Wdiv>0);  % find divergent points in space and time
	 fpd = intersect(fd,fp); % find cells where >15% ice and divergent
         fdf = intersect(fd,ffrz); % freezing and divergent points
         fpf = intersect(fp,ffrz); % points where ice >15% and freezing
        % area-weighted growth rate only over cells where ice is growing
         GR_20C(n,k) = sum((WI2O(ffrz).*rep_cellareaxx(ffrz))/sum(rep_cellareaxx(ffrz)));
        % total ice production (m^3). Need factor 30 to convert from mean daily to mean monthly 
	 iceprod_20C(n,k) = sum(WI2O(ffrz).*rep_cellareaxx(ffrz))*30;
	 % total positive area diverged each winter in m^2
	 pos_areadiv_20C(n,k) = sum((Wdiv(fd)*(30/100)).*rep_cellareaxx(fd)); % factor (30/100) converts %/day to /month
	 % net area diverged each winter in m^2
	 net_areadiv_20C(n,k) = sum((Wdiv(fnd)*(30/100)).*rep_cellareaxx(fnd)); % factor (30/100) converts %/day to /month
         % snow - count over >15% ice cells
         hs_20C(n,k) = sum((Wsnow(fp).*rep_cellareaxx(fp))/sum(rep_cellareaxx(fp)));
	 % September sea ice area
	 SepSIA_20C(n,k) = squeeze(nansum(nansum(SepSIC.*cellareaxx)));
         % temperature - count over all ocean cells
	 SAT_20C(n,k) = sum((WSAT(fn).*rep_cellareaxx(fn))/sum(rep_cellareaxx(fn)));
         % freezing area days
         M1 = ones(size(WI2O))*30;
         frz_area_days_20C(n,k) = sum(M1(ffrz).*rep_cellareaxx(ffrz));
         % heat fluxes - count and normalise over days where ice present at >15%
         Fw_20C(n,k) = sum((WFw(ffrz).*rep_cellareaxx(ffrz))/sum(rep_cellareaxx(ffrz)));


        for n = 1:73 % no. of years mins 2 incomplete winter seasons
                 Wdiv = NaN*zeros(size(div,1),size(div,2),7);
                 WI2O = NaN*zeros(size(I2O,1),size(I2O,2),7);
		 WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
                 WSepSIC = NaN*zeros(size(SIC,1),size(SIC,2));
                 Wsnow = NaN*zeros(size(snow,1),size(snow,2),7);
                 WSAT = NaN*zeros(size(SAT,1),size(SAT,2),7);
                 WFw = NaN*zeros(size(SIC,1),size(SIC,2),7);
         % take winter months
         Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);
         WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16).*0.01;% assuming the I2O has ice conc built into it, converting cm /day --> m/day
	 WSIC = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+16);
         SepSIC = msk_SIC(:,:,12*(n-1)+9);
	 Wsnow = msk_snow(:,:,12*(n-1)+10:12*(n-1)+16);
         WSAT = msk_SAT(:,:,12*(n-1)+10:12*(n-1)+16);
         WFw = msk_Fw(:,:,12*(n-1)+10:12*(n-1)+16);
         
	 fn = find(~isnan(WI2O)); % find non NaN values; this finds values over the ocean - however it looks like SAT is also NaN over land as from cice
         fp = find(WSIC>0.15); % find cells with more than 15% sea ice
	 fnd = find(~isnan(Wdiv)); % find non NaN values of divergence
         ffrz = find(WI2O<0); % find freezing points in time and space
         fd = find(Wdiv>0);  % find divergent points in space and time
	 fpd = intersect(fd,fp); % find cells where >15% ice and divergent
         fdf = intersect(fd,ffrz); % freezing and divergent points
         fpf = intersect(fp,ffrz); % points where ice >15% and freezing
        % area-weighted growth rate only over cells where ice is growing
         GR_RCP85(n,k) = sum((WI2O(ffrz).*rep_cellareaxx(ffrz))/sum(rep_cellareaxx(ffrz)));
        % total ice production (m^3). Need factor 30 to convert from mean daily to mean monthly. 
	 iceprod_RCP85(n,k) = sum(WI2O(ffrz).*rep_cellareaxx(ffrz))*30;
	 % total positive area diverged each winter in m^2
	 pos_areadiv_RCP85(n,k) = sum((Wdiv(fd)*(30/100)).*rep_cellareaxx(fd)); % factor (30/100) converts %/day to /month
	 % net area diverged each winter in m^2
	 net_areadiv_RCP85(n,k) = sum((Wdiv(fnd)*(30/100)).*rep_cellareaxx(fnd)); % factor (30/100) converts %/day to /month
         % snow - count over >15% ice cells
         hs_RCP85(n,k) = sum((Wsnow(fp).*rep_cellareaxx(fp))/sum(rep_cellareaxx(fp)));
	 % September sea ice area
	 SepSIA_RCP85(n,k) = squeeze(nansum(nansum(SepSIC.*cellareaxx)));
         % temperature - count over all ocean cells
         SAT_RCP85(n,k) = sum((WSAT(fn).*rep_cellareaxx(fn))/sum(rep_cellareaxx(fn)));
         % freezing area days
         M1 = ones(size(WI2O))*30;
         frz_area_days_RCP85(n,k) = sum(M1(ffrz).*rep_cellareaxx(ffrz));
	 % heat fluxes - count and normalise over days where ice present at >15%
	 Fw_RCP85(n,k) = sum((WFw(ffrz).*rep_cellareaxx(ffrz))/sum(rep_cellareaxx(ffrz)));
	 end     
 fprintf('finished RCP85 ensemble # %d\n',k)
end

save /home/ocean_personal_data/samc/CESM/variables/CESM_iceprod_vars_36to40.mat GR_20C GR_RCP85 iceprod_20C iceprod_RCP85 net_areadiv_20C net_areadiv_RCP85 pos_areadiv_20C pos_areadiv_RCP85 hs_20C hs_RCP85 SepSIA_20C SepSIA_RCP85 SAT_20C SAT_RCP85 frz_area_days_20C frz_area_days_RCP85 Fw_20C Fw_RCP85

