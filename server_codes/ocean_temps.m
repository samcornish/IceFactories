% Script to retrieve 10 m ocean temperature
clear
close all
addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/TEMP/
addpath /home/ocean_shared_data2/CESM_arctic/
addpath /home/ocean_personal_data/samc/CESM/variables/
addpath /home/ocean_shared_data2/CESM_arctic/fresh/

load KLmask_pop.mat
load fnames40.mat

cellarea = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.TEMP.192001-200512.nc','TAREA');
rep_cellarea = repmat(cellarea,[1,1,8]); % repeat the matrix 8 times in a third dimension in order to allow multiplying with winter arrays pointwise
z = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.TEMP.192001-200512.nc','z_t')/100;	% depth to midpoint of each cell in m
dz = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.TEMP.192001-200512.nc','dz')/100; % thickness of each cell
cellvol = repmat(cellarea,[1,1,5])*10;
rep_cellvol = repmat(cellvol,[1,1,1,8]);

monthly_T10m_20C = zeros(84,8,32); % initialise, start from Septerm *different to other monthly records*
monthly_T50m_20C= zeros(84,8,32); % initialise, start from Septerm *different to other monthly records*
 
for k = 1:35	% total number of ensembles that we will use - member 001 is missing. If fix this, need to correct indices in script below.
% 20th CENTURY
if k == 1 
	T10m = squeeze(ncread(strcat(fnames_20C{k},'.pop.h.TEMP.185001-200512.nc'),'TEMP',[1,1,1,841],[320,49,1,1032]));
        T50m = squeeze(ncread(strcat(fnames_20C{k},'.pop.h.TEMP.185001-200512.nc'),'TEMP',[1,1,1,841],[320,49,5,1032]));
else
	T10m = squeeze(ncread(strcat(fnames_20C{k},'.pop.h.TEMP.192001-200512.nc'),'TEMP',[1,1,1,1],[320,49,1,1032]));
    	T50m = squeeze(ncread(strcat(fnames_20C{k},'.pop.h.TEMP.192001-200512.nc'),'TEMP',[1,1,1,1],[320,49,5,1032]));
end
msk_T10m = zeros(size(T10m));
msk_T50m = zeros(size(T50m));
for i = 1:size(T10m,3)
msk_T10m(:,:,i) = T10m(:,:,i).*mask;
msk_T50m(:,:,:,i) = T50m(:,:,:,i).*repmat(mask,[1,1,5]);
end

for n = 1:85 % no of years - 1 incomplete winter season
    WT10m = NaN*zeros(size(T10m,1),size(T10m,2),8);
    WT50m = NaN*zeros(size(T10m,1),size(T10m,2),size(T50m,3),8);
    % take winter month but start in Seps
    WT10m = msk_T10m(:,:,12*(n-1)+9:12*(n-1)+16);
    WT50m = msk_T50m(:,:,:,12*(n-1)+9:12*(n-1)+16);
    if n == 1
    fn = find(~isnan(WT10m)); % find non NaN values; this finds values over the ocean
	M = zeros(size(WT10m));
	M(fn) = 1; 
    fn50 = find(~isnan(WT50m)); %find non NaN values in top 50m
    	M50 = zeros(size(WT50m));
	M50(fn50) = 1; % depth of every cell in top 100 m is 10 m 
    	vol_sum50 = squeeze(nansum(nansum(nansum(M50.*rep_cellvol))));
	area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));	% sum of all areas where T10m is recorded in mask
    end
    	monthly_T10m_20C(n,:,k) = squeeze(nansum(nansum(WT10m.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions and give us the area-weighted 10 m temperature
	monthly_T50m_20C(n,:,k) = squeeze(nansum(nansum(nansum(WT50m.*rep_cellvol))))./vol_sum50;
end   

% now save the first part of the crossover year
n = 86;
    WT10m = NaN*zeros(size(T10m,1),size(T10m,2),8);
    WT50m = NaN*zeros(size(T10m,1),size(T10m,2),size(T50m,3),8);
    % take winter month but start in Seps
    WT10m(:,:,1:4) = msk_T10m(:,:,12*(n-1)+9:12*(n-1)+12);
    WT50m(:,:,:,1:4) = msk_T50m(:,:,:,12*(n-1)+9:12*(n-1)+12);


% RCP8.5 
if k >= 34
	T10m = squeeze(ncread(strcat(fnames_RCP85{k},'.pop.h.TEMP.200601-210012.nc'),'TEMP',[1,1,1,1],[320,49,1,900]));
        T50m = squeeze(ncread(strcat(fnames_RCP85{k},'.pop.h.TEMP.200601-210012.nc'),'TEMP',[1,1,1,1],[320,49,5,900]));

else
	T10m = squeeze(ncread(strcat(fnames_RCP85{k},'.pop.h.TEMP.200601-208012.nc'),'TEMP',[1,1,1,1],[320,49,1,900]));
	T50m = squeeze(ncread(strcat(fnames_RCP85{k},'.pop.h.TEMP.200601-208012.nc'),'TEMP',[1,1,1,1],[320,49,5,900]));
end
msk_T10m = zeros(size(T10m));
msk_T50m = zeros(size(T50m));

for i = 1:size(T10m,3)
msk_T10m(:,:,i) = T10m(:,:,i).*mask;
msk_T50m(:,:,:,i) = T50m(:,:,:,i).*repmat(mask,[1,1,5]);
end

% first of all, process the first year - note that this will be associated with 20C as the last year.
    % take winter months
    WT10m(:,:,5:8) = msk_T10m(:,:,1:4);
    WT50m(:,:,:,5:8) = msk_T50m(:,:,:,1:4);
    if n == 1
    fn = find(~isnan(WT10m)); % find non NaN values; this finds values over the ocean
        M = zeros(size(WT10m));
        M(fn) = 1;
    fn50 = find(~isnan(WT50m)); %find non NaN values in top 50m
        M50 = zeros(size(WT50m));
        M50(fn50) = 1; % depth of every cell in top 100 m is 10 m 
        vol_sum50 = squeeze(nansum(nansum(nansum(M50.*rep_cellvol))));
        area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));    % sum of all areas where T10m is recorded in mask
    end
        monthly_T10m_20C(n,:,k) = squeeze(nansum(nansum(WT10m.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions and give us the area-weighted 10 m temperature
        monthly_T50m_20C(n,:,k) = squeeze(nansum(nansum(nansum(WT50m.*rep_cellvol))))./vol_sum50;

fprintf('finished 20C ensemble # %d\n',k)


for n = 1:73 % no of years - 2 incomplete winter seasons
    WT10m = NaN*zeros(size(T10m,1),size(T10m,2),8);
    WT50m = NaN*zeros(size(T10m,1),size(T10m,2),size(T50m,3),8);
    % take winter month but start in Seps
    WT10m = msk_T10m(:,:,12*(n-1)+9:12*(n-1)+16);
    WT50m = msk_T50m(:,:,:,12*(n-1)+9:12*(n-1)+16);
    if n==1
    fn = find(~isnan(WT10m)); % find non NaN values; this finds values over the ocean
	M = zeros(size(WT10m));
	M(fn) = 1; 
    fn50 = find(~isnan(WT50m)); %find non NaN values in top 50m
        M50 = zeros(size(WT50m));
        M50(fn50) = 1; % depth of every cell in top 100 m is 10 m 
        vol_sum50 = squeeze(nansum(nansum(nansum(M50.*rep_cellvol))));
	area_sum = squeeze(nansum(nansum(M.*rep_cellarea)));	% sum of all areas where T10m is recorded in mask
    end
    	monthly_T10m_RCP85(n,:,k) = squeeze(nansum(nansum(WT10m.*rep_cellarea)))./area_sum; % this should mean through the first then second dimensions and give us the area-weighted 10 m temperature
	monthly_T50m_RCP85(n,:,k) = squeeze(nansum(nansum(nansum(WT50m.*rep_cellvol))))./vol_sum50;
end   

fprintf('finished RCP85 ensemble # %d\n',k)
end

save /home/ocean_personal_data/samc/CESM/variables/CESM_T10m_20C_RCP85.mat monthly_T10m_20C monthly_T10m_RCP85 monthly_T50m_20C monthly_T50m_RCP85  

