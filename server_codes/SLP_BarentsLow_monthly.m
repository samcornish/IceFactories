% script to take monthly values of SLP at the point where R^2 with div is maximised (a point in the Barents Sea)
% I will then use these SLP values to correlate against SIC, to see when a correlation emerges (and in which months)

clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/PSL/
addpath /home/ocean_personal_data/samc/CESM/variables/

load fnames_div.mat
load slp_mpd_regression.mat R2

[val ind] = max(R2(:));
[lowi lowj] = ind2sub([size(R2,1),size(R2,2)],ind);
lowi = lowi+144;
lat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.PSL.192001-200512.nc','lat');
lon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cam.h0.PSL.192001-200512.nc','lon');

for k = 1:32	% total number of ensembles that we will use
% 20th CENTURY
	SLP = squeeze(ncread(strcat(fnames_div_20C{k},'.cam.h0.PSL.192001-200512.nc'),'PSL',[lowi,lowj,1],[1,1,1032]));

for n = 1:84 % no of years - 2 incomplete winter seasons
    Wslp = zeros(1,7);
    % take winter months
    Wslp = SLP(12*(n-1)+10:12*(n-1)+16);
    monthly_lowP_20C(n,:,k) = Wslp;
end
    fprintf('finished 20C ensemble # %d\n',k)

    clearvars SLP 

% RCP8.5 
	SLP = squeeze(ncread(strcat(fnames_div_RCP85{k},'.cam.h0.PSL.200601-208012.nc'),'PSL',[lowi,lowj,1],[1,1,900]));

for n = 1:73 % no of years - 2 incomplete winter seasons
    Wslp = zeros(1,7);
    % take winter months
    Wslp = SLP(12*(n-1)+10:12*(n-1)+16);
    monthly_lowP_RCP85(n,:,k) = Wslp;
end

fprintf('finished RCP8.5 ensemble # %d\n',k)
end
save /home/ocean_personal_data/samc/CESM/variables/SLP_BarentsLow_monthly.mat monthly_lowP_20C monthly_lowP_RCP85 lat lon lowi lowj val 
