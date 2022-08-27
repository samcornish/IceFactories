% script to correlate div and SIA with one another at each month, through the ensemble

clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/divu/
addpath /home/ocean_shared_data2/CESM_arctic/aice/
addpath /home/ocean_personal_data/samc/CESM/variables/

load KLmask.mat
load fnames_div.mat

dyt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dyt');
dxt = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','dxt');
cellarea = dyt.*dxt;
monthly_SIA_20C = zeros(84,7,32); % initialise
monthly_SIA_RCP85 = zeros(73,7,32); % initialise

for k = 1:32	% total number of ensembles that we will use
% 20th CENTURY
	div = ncread(strcat(fnames_div_20C{k},'.cice.h.divu_nh.192001-200512.nc'),'divu');
	SIC = ncread(strcat(fnames_div_20C{k},'.cice.h.aice_nh.192001-200512.nc'),'aice')/100;
msk_SIC = zeros(size(SIC));
msk_div = zeros(size(div));
for i = 1:size(SIC,3)
msk_SIC(:,:,i) = SIC(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end

for n = 1:84 % no of years - 2 incomplete winter seasons
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
    Wdiv = NaN*zeros(size(div,1),size(div,2),7);

    % take winter months
    WSIC = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+16);
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);
   % take spatial means for each month
   smWSIC = squeeze(nanmean(nanmean(WSIC)));
   smWdiv = squeeze(nanmean(nanmean(Wdiv))); 
   monthly_SIC_20C(n,:,k) = smWSIC;
   monthly_div_20C(n,:,k) = smWdiv;
end   

fprintf('finished 20C ensemble # %d\n',k)

% RCP8.5 

	SIC = ncread(strcat(fnames_div_RCP85{k},'.cice.h.aice_nh.200601-208012.nc'),'aice')/100;
	div = ncread(strcat(fnames_div_RCP85{k},'.cice.h.divu_nh.200601-208012.nc'),'divu');

msk_SIC = zeros(size(SIC));
msk_div = zeros(size(div));
for i = 1:size(SIC,3)
msk_SIC(:,:,i) = SIC(:,:,i).*mask;
msk_div(:,:,i) = div(:,:,i).*mask;
end
for n = 1:73 % no of years - 2 incomplete winter seasons
    WSIC = NaN*zeros(size(SIC,1),size(SIC,2),7);
    Wdiv = NaN*zeros(size(div,1),size(div,2),7);

    % take winter months
    WSIC = msk_SIC(:,:,12*(n-1)+10:12*(n-1)+16);
    Wdiv = msk_div(:,:,12*(n-1)+10:12*(n-1)+16);
   % take spatial means for each month
   smWSIC = squeeze(nanmean(nanmean(WSIC)));
   smWdiv = squeeze(nanmean(nanmean(Wdiv))); 
   monthly_SIC_RCP85(n,:,k) = smWSIC;
   monthly_div_RCP85(n,:,k) = smWdiv;
end   

fprintf('finished RCP85 ensemble # %d\n',k)
end
%% loop through years and take correlations through the ensemble

load SLP_BarentsLow_monthly.mat monthly_lowP_20C monthly_lowP_RCP85

for n = 1:84
	for j = 1:7
		[r p] = corrcoef(monthly_SIC_20C(n,j,:),monthly_div_20C(n,j,:));
		R_monthlySICdiv_20C(n,j) = r(1,2);
		P_monthlySICdiv_20C(n,j) = p(1,2);
		[r p] = corrcoef(monthly_SIC_20C(n,j,:),monthly_lowP_20C(n,j,:));
		R_monthlySIClowP_20C(n,j) = r(1,2);
		P_monthlySIClowP_20C(n,j) = p(1,2);
	end
end


for n = 1:73
	for j = 1:7
		[r p] = corrcoef(monthly_SIC_RCP85(n,j,:),monthly_div_RCP85(n,j,:));
		R_monthlySICdiv_RCP85(n,j) = r(1,2);
		P_monthlySICdiv_RCP85(n,j) = p(1,2);
		[r p] = corrcoef(monthly_SIC_RCP85(n,j,:),monthly_lowP_RCP85(n,j,:));
		R_monthlySIClowP_RCP85(n,j) = r(1,2);
		P_monthlySIClowP_RCP85(n,j) = p(1,2);
	end
end

save /home/ocean_personal_data/samc/CESM/variables/CESM_SIC_div_monthly.mat monthly_SIC_20C monthly_div_20C monthly_SIC_RCP85 monthly_div_RCP85 monthly_lowP_20C monthly_lowP_RCP85 R_monthlySICdiv_20C P_monthlySICdiv_20C R_monthlySICdiv_RCP85 P_monthlySICdiv_RCP85 R_monthlySIClowP_20C P_monthlySIClowP_20C R_monthlySIClowP_RCP85 P_monthlySIClowP_RCP85 
