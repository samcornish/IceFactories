% Script to retrieve winter mean air temperatures, monthly air temperatures, and perform multi regression of monthly air T against total ice production 
clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/Tair/
addpath /home/ocean_shared_data2/CESM_arctic/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/
addpath /home/ocean_shared_data2/CESM_arctic/divu/

load KLmask.mat
load fnames40.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.Tair_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.Tair_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
% Now using 7 months starting Oct - saving variable under different filename
monthly_SAT_20C = zeros(84,7,32); % initialise
monthly_SAT_RCP85 = zeros(73,7,32);
rep_cellarea = repmat(cellarea,[1,1,7]); % repeat the matrix 7 times in a third dimension in order to allow multiplying with winter arrays pointwise
pcFRZ = 50/100; %percent freezing - threshold for taking SAT and onset timing

for k = 1:35	% total number of ensembles that we will use
% 20th CENTURY
if k == 1
	SAT = ncread(strcat(fnames_20C{k},'.cice.h.Tair_nh.185001-200512.nc'),'Tair',[1,1,841],[291,48,1032]);
        I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.185001-200512.nc'),'fresh',[1,1,841],[291,48,1032]);
        div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.185001-200512.nc'),'divu',[1,1,841],[291,48,1032]);
else
	SAT = ncread(strcat(fnames_20C{k},'.cice.h.Tair_nh.192001-200512.nc'),'Tair');
	I2O = ncread(strcat(fnames_20C{k},'.cice.h.fresh_nh.192001-200512.nc'),'fresh');
	div = ncread(strcat(fnames_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
end
msk_SAT = zeros(size(SAT));
msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
for i = 1:size(SAT,3)
msk_SAT(:,:,i) = SAT(:,:,i).*mask;
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end

for n = 1:85 % no of years - 2 incomplete winter seasons
    WSAT = NaN*zeros(size(SAT,1),size(SAT,2),7);
    WI2O = zeros(size(I2O,1),size(I2O,2),7);
    Wdiv = zeros(size(div,1),size(div,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WSAT = msk_SAT(:,:,12*(n-1)+10:12*(n-1)+16);
    fn = find(~isnan(WI2O)); % find non NaN values; this finds values over the ocean - however it looks like SAT is also NaN over land as from cice
    ffrz = find(WI2O<0); % find freezing points in time and space
    fnfrz = find(WI2O>=0); % find non-freezing points in time and space
    fqd = find(Wdiv<0.5); fqc = find(Wdiv>-0.5); fq = intersect(fqd,fqc); % find points that are quiescent with respect to divergence
    SAT_20C(n,k) = sum((WSAT(fn).*rep_cellarea(fn))/sum(rep_cellarea(fn))); % find the area-weighted mean SAT 
    SATq_20C(n,k) = sum((WSAT(fq).*rep_cellarea(fq))/sum(rep_cellarea(fq))); % find the area-weighted mean SAT
    SAT_frz_20C(n,k) = sum((WSAT(ffrz).*rep_cellarea(ffrz))/sum(rep_cellarea(ffrz))); % find the area-weighted mean SAT on freezing points
    SAT_nfrz_20C(n,k) = sum((WSAT(fnfrz).*rep_cellarea(fnfrz))/sum(rep_cellarea(fnfrz))); % find the area-weighted mean SAT where not freezing	
	M = zeros(size(WSAT));
	M(fn) = 1; 
	area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));	% sum of all areas where SAT is recorded in mask
    monthly_SAT_20C(n,:,k) = squeeze(nansum(nansum(WSAT.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    WSATm = zeros(size(WSAT));
    for m = 1:5	% go through each month up to mid-winter: don't want April showing up
	    fnm = find(~isnan(WI2O(:,:,m))); % find non NaN values that month;
	    mfrz = find(WI2O(:,:,m)<0); % find freezing points that month
	if length(mfrz)<pcFRZ*length(fnm) % if there are fewer than 90% points freezing, continue to populate this matrix
		WSATm(:,:,m) = WSAT(:,:,m);
		onset_month_20C(n,k) = m;
	elseif m  == 1	% always populate using the first month - may be wise to start in Sep for this
		WSATm(:,:,m) = WSAT(:,:,m);
		onset_month_20C(n,k) = m;
	else	% if not first month and more than the threshold number of points are freezing, do this:
		sumWSATm = sum(WSATm,3)/m-1;	% sum through and divide by number of months used (mean)	
		SAT_mb4frz_20C(n,k) =sum((sumWSATm(fnm).*cellarea(fnm))/sum(cellarea(fnm))); % find the area-weighted mean SAT in the months before freezing threshold
		break % break the loop in m if there are more than the threshold points freezing after the first month 
	end   
	end
end

fprintf('finished 20C ensemble # %d\n',k)


% now save the first part of the crossover year
n = 86;
    WSAT = NaN*zeros(size(SAT,1),size(SAT,2),7);
    WI2O = zeros(size(I2O,1),size(I2O,2),7);
    Wdiv = zeros(size(div,1),size(div,2),7);
    % take winter months
    Wdiv(:,:,1:3) = msk_div(:,:,12*(n-1)+10:12*(n-1)+12);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+12);
    WI2O(:,:,1:3) = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+12);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+12);
    WSAT(:,:,1:3) = msk_SAT(:,:,12*(n-1)+10:12*(n-1)+12);

% RCP8.5 
if k >= 34
 	I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-210012.nc'),'fresh',[1,1,1],[291,48,900]);
        SAT = ncread(strcat(fnames_RCP85{k},'.cice.h.Tair_nh.200601-210012.nc'),'Tair',[1,1,1],[291,48,900]);
        div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-210012.nc'),'divu',[1,1,1],[291,48,900]);
else
        I2O = ncread(strcat(fnames_RCP85{k},'.cice.h.fresh_nh.200601-208012.nc'),'fresh');
        SAT = ncread(strcat(fnames_RCP85{k},'.cice.h.Tair_nh.200601-208012.nc'),'Tair');
	div = ncread(strcat(fnames_RCP85{k},'.cice.h.divu_nh.200601-208012.nc'),'divu');
end
msk_SAT = zeros(size(SAT));
msk_I2O = zeros(size(I2O));
msk_div = zeros(size(div));
for i = 1:size(SAT,3)
msk_SAT(:,:,i) = SAT(:,:,i).*mask;
msk_I2O(:,:,i) = I2O(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end

% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months

    Wdiv(:,:,4:7) = msk_div(:,:,1:4);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+12);
    WI2O(:,:,4:7) = msk_I2O(:,:,1:4);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+12);
    WSAT(:,:,4:7) = msk_SAT(:,:,1:4);
    fn = find(~isnan(WI2O)); % find non NaN values; this finds values over the ocean - however it looks like SAT is also NaN over land as from cice
    ffrz = find(WI2O<0); % find freezing points in time and space
    fnfrz = find(WI2O>=0); % find non-freezing points in time and space
    fqd = find(Wdiv<0.5); fqc = find(Wdiv>-0.5); fq = intersect(fqd,fqc); % find points that are quiescent with respect to divergence
    SAT_20C(n,k) = sum((WSAT(fn).*rep_cellarea(fn))/sum(rep_cellarea(fn))); % find the area-weighted mean SAT 
    SATq_20C(n,k) = sum((WSAT(fq).*rep_cellarea(fq))/sum(rep_cellarea(fq))); % find the area-weighted mean SAT
    SAT_frz_20C(n,k) = sum((WSAT(ffrz).*rep_cellarea(ffrz))/sum(rep_cellarea(ffrz))); % find the area-weighted mean SAT on freezing points
    SAT_nfrz_20C(n,k) = sum((WSAT(fnfrz).*rep_cellarea(fnfrz))/sum(rep_cellarea(fnfrz))); % find the area-weighted mean SAT where not freezing
        M = zeros(size(WSAT));
        M(fn) = 1;
        area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));    % sum of all areas where SAT is recorded in mask
    monthly_SAT_20C(n,:,k) = squeeze(nansum(nansum(WSAT.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    WSATm = zeros(size(WSAT));
    for m = 1:5 % go through each month up to mid-winter: don't want April showing up
            fnm = find(~isnan(WI2O(:,:,m))); % find non NaN values that month;
            mfrz = find(WI2O(:,:,m)<0); % find freezing points that month
        if length(mfrz)<pcFRZ*length(fnm) % if there are fewer than 90% points freezing, continue to populate this matrix
                WSATm(:,:,m) = WSAT(:,:,m);
                onset_month_20C(n,k) = m;
        elseif m  == 1  % always populate using the first month - may be wise to start in Sep for this
                WSATm(:,:,m) = WSAT(:,:,m);
                onset_month_20C(n,k) = m;
        else    % if not first month and more than the threshold number of points are freezing, do this:
                sumWSATm = sum(WSATm,3)/m-1;    % sum through and divide by number of months used (mean)
                SAT_mb4frz_20C(n,k) =sum((sumWSATm(fnm).*cellarea(fnm))/sum(cellarea(fnm))); % find the area-weighted mean SAT in the months before freezing threshold
                break % break the loop in m if there are more than the threshold points freezing after the first month 
        end
        end


for n = 1:73 % no of years - 2 incomplete winter seasons
    WSAT = NaN*zeros(size(SAT,1),size(SAT,2),7);
    WI2O = zeros(size(I2O,1),size(I2O,2),7);
    Wdiv = zeros(size(div,1),size(div,2),7);
    % take winter months
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WI2O = msk_I2O(:,:,12*(n-1)+10:12*(n-1)+16);%.*SIconc_rg(:,:,12*(n-1)+10:12*(n-1)+16);
    WSAT = msk_SAT(:,:,12*(n-1)+10:12*(n-1)+16);
    fn = find(~isnan(WI2O)); % find non NaN values; this finds values over the ocean - however it looks like SAT is also NaN over land as from cice
    ffrz = find(WI2O<0); % find freezing points in time and space
    fnfrz = find(WI2O>=0); % find non-freezing points in time and space
    fqd = find(Wdiv<0.5); fqc = find(Wdiv>-0.5); fq = intersect(fqd,fqc); % find points that are quiescent with respect to divergence
    SAT_RCP85(n,k) = sum((WSAT(fn).*rep_cellarea(fn))/sum(rep_cellarea(fn))); % find the area-weighted mean SAT 
    SATq_RCP85(n,k) = sum((WSAT(fq).*rep_cellarea(fq))/sum(rep_cellarea(fq))); % find the area-weighted mean SAT 
    SAT_frz_RCP85(n,k) = sum((WSAT(ffrz).*rep_cellarea(ffrz))/sum(rep_cellarea(ffrz))); % find the area-weighted mean SAT on freezing points
    SAT_nfrz_RCP85(n,k) = sum((WSAT(fnfrz).*rep_cellarea(fnfrz))/sum(rep_cellarea(fnfrz))); % find the area-weighted mean SAT where not freezing	
	M = zeros(size(WSAT));
	M(fn) = 1; 
	area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));	% sum of all areas where SAT is recorded in mask
    monthly_SAT_RCP85(n,:,k) = squeeze(nansum(nansum(WSAT.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions
    WSATm = zeros(size(WSAT));
    for m = 1:5	% go through each month up to mid-winter: don't want April showing up
	    fnm = find(~isnan(WI2O(:,:,m))); % find non NaN values that month;
	    mfrz = find(WI2O(:,:,m)<0); % find freezing points that month
	if length(mfrz)<pcFRZ*length(fnm) % if there are fewer than 90% points freezing, continue to populate this matrix
		WSATm(:,:,m) = WSAT(:,:,m);
		onset_month_RCP85(n,k) = m;
	elseif m  == 1	% always populate using the first month - may be wise to start in Sep for this
		WSATm(:,:,m) = WSAT(:,:,m);
		onset_month_RCP85(n,k) = m;
	else	% if not first month and more than the threshold number of points are freezing, do this:
		sumWSATm = sum(WSATm,3)/m-1;	% sum through and divide by number of months used (mean)	
		SAT_mb4frz_RCP85(n,k) =sum((sumWSATm(fnm).*cellarea(fnm))/sum(cellarea(fnm))); % find the area-weighted mean SAT in the months before freezing threshold
		break % break the loop in m if there are more than the threshold points freezing after the first month 
	end   
	end
end   

fprintf('finished RCP8.5 ensemble # %d\n',k)
end

%% now undertake the monthly regressions
%load CESM_div_I2O_reg_20C_RCP85.mat iceprod_20C 
%for i = 1:7
%M_SAT = reshape(monthly_SAT_20C(1:60,i,:),60*32,1);	% reshape monthly SAT
%X = [ones(60*32,1),M_SAT];	% append a column of ones to the beginning
%Y = reshape(iceprod_20C(1:60,:),60*32,1)/1e11; % divide by 1e11 to get X and Y on same order of magnitude
%
%[b,bint,r,rint,stats] = regress(Y,X);
%monthly_simplereg_SAT_beta(i,:) = b*1e11;
%monthly_simpleR2(i) = stats(1);
%end

save /home/ocean_personal_data/samc/CESM/variables/CESM_Tair_20C_RCP85.mat SAT_20C SATq_20C SAT_RCP85 SATq_RCP85 monthly_SAT_20C monthly_SAT_RCP85 SAT_frz_20C SAT_nfrz_20C SAT_frz_RCP85 SAT_nfrz_RCP85 onset_month_20C onset_month_RCP85 SAT_mb4frz_20C SAT_mb4frz_RCP85 
