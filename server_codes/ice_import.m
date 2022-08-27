% script to calculate ice volume export from the Kara Laptev Seas

clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/uvel/
addpath /home/ocean_shared_data2/CESM_arctic/vvel/
addpath /home/ocean_shared_data2/CESM_arctic/hi/
addpath /home/ocean_shared_data2/CESM_arctic/aice/
addpath /home/ocean_personal_data/samc/CESM/variables/

load KLmask.mat
load fnames_div.mat

[fix fjx] = find(xmask);
[fiy fjy] = find(ymask);	% locations of export mask (on edge of region)
xmask_m = xmask.*dyu;
ymask_m = ymask.*dxu;
for k = 1:32    % total number of ensembles that we will use
% 20th CENTURY
        v = ncread(strcat(fnames_div_20C{k},'.cice.h.vvel_nh.192001-200512.nc'),'vvel');
        u = ncread(strcat(fnames_div_20C{k},'.cice.h.uvel_nh.192001-200512.nc'),'uvel');
        hi = ncread(strcat(fnames_div_20C{k},'.cice.h.hi_nh.192001-200512.nc'),'hi');
        SIC = ncread(strcat(fnames_div_20C{k},'.cice.h.aice_nh.192001-200512.nc'),'aice')/100;
Whi_int = zeros(size(xmask));
WSIC_int = zeros(size(xmask));
for n = 1:84
   Wv = zeros(size(v,1),size(v,2),7);
   Wu = zeros(size(v,1),size(v,2),7);
   Whi = zeros(size(hi,1),size(hi,2),7);
   WSIC = zeros(size(SIC,1),size(SIC,2),7);
    % take winter months
    Wv = v(:,:,12*(n-1)+10:12*(n-1)+16);
    Wu = u(:,:,12*(n-1)+10:12*(n-1)+16);
    Whi = hi(:,:,12*(n-1)+10:12*(n-1)+16);
    WSIC = SIC(:,:,12*(n-1)+10:12*(n-1)+16);
    for t = 1:7
    for i = 1:length(fix)
        Whi_int = nanmean([Whi(fix(i)+4,fjx(i)+1,t);Whi(fix(i)+4,fjx(i)+2,t);Whi(fix(i)+5,fjx(i)+1,t);Whi(fix(i)+5,fjx(i)+2,t)]);
        WSIC_int = nanmean([WSIC(fix(i)+4,fjx(i)+1,t);WSIC(fix(i)+4,fjx(i)+2,t);WSIC(fix(i)+5,fjx(i)+1,t);WSIC(fix(i)+5,fjx(i)+2,t)]);
        volumeTx(i) = Whi_int*xmask_m(fix(i),fjx(i))*Wu(fix(i),fjx(i),t)*30*86400/100;
        areaTx(i) = WSIC_int*xmask_m(fix(i),fjx(i))*Wu(fix(i),fjx(i),t)*30*86400/100;
    end
    for i = 1:length(fiy)
        Whi_int = nanmean([Whi(fiy(i)+4,fjy(i)+1,t);Whi(fiy(i)+4,fjy(i)+2,t);Whi(fiy(i)+5,fjy(i)+1,t);Whi(fiy(i)+5,fjy(i)+2,t)]);
        WSIC_int = nanmean([WSIC(fiy(i)+4,fjy(i)+1,t);WSIC(fiy(i)+4,fjy(i)+2,t);WSIC(fiy(i)+5,fjy(i)+1,t);WSIC(fiy(i)+5,fjy(i)+2,t)]);
        volumeTy(i) = Whi_int*ymask_m(fiy(i),fjy(i))*Wv(fiy(i),fjy(i),t)*30*86400/100;
        areaTy(i) = WSIC_int*ymask_m(fiy(i),fjy(i))*Wv(fiy(i),fjy(i),t)*30*86400/100;
    end
    Vimport_monthly_20C(n,t,k) = nansum(volumeTx) + nansum(volumeTy);
    Aimport_monthly_20C(n,t,k) = nansum(areaTx) + nansum(areaTy);
	clearvars volumeTx volumeTy areaTx areaTy
    end
end
    fprintf('finished 20C ensemble # %d\n',k)

% RCP85

 	v = ncread(strcat(fnames_div_RCP85{k},'.cice.h.vvel_nh.200601-208012.nc'),'vvel');
        u = ncread(strcat(fnames_div_RCP85{k},'.cice.h.uvel_nh.200601-208012.nc'),'uvel');
        hi = ncread(strcat(fnames_div_RCP85{k},'.cice.h.hi_nh.200601-208012.nc'),'hi');
        SIC = ncread(strcat(fnames_div_RCP85{k},'.cice.h.aice_nh.200601-208012.nc'),'aice')/100;
Whi_int = zeros(size(xmask));
WSIC_int = zeros(size(xmask));
for n = 1:73
   Wv = zeros(size(v,1),size(v,2),7);
   Wu = zeros(size(v,1),size(v,2),7);
   Whi = zeros(size(hi,1),size(hi,2),7);
   WSIC = zeros(size(SIC,1),size(SIC,2),7);
    % take winter months
    Wv = v(:,:,12*(n-1)+10:12*(n-1)+16);
    Wu = u(:,:,12*(n-1)+10:12*(n-1)+16);
    Whi = hi(:,:,12*(n-1)+10:12*(n-1)+16);
    WSIC = SIC(:,:,12*(n-1)+10:12*(n-1)+16);
    for t = 1:7
    for i = 1:length(fix)
        Whi_int = nanmean([Whi(fix(i)+4,fjx(i)+1,t);Whi(fix(i)+4,fjx(i)+2,t);Whi(fix(i)+5,fjx(i)+1,t);Whi(fix(i)+5,fjx(i)+2,t)]);
        WSIC_int = nanmean([WSIC(fix(i)+4,fjx(i)+1,t);WSIC(fix(i)+4,fjx(i)+2,t);WSIC(fix(i)+5,fjx(i)+1,t);WSIC(fix(i)+5,fjx(i)+2,t)]);
        volumeTx(i) = Whi_int*xmask_m(fix(i),fjx(i))*Wu(fix(i),fjx(i),t)*30*86400/100;
        areaTx(i) = WSIC_int*xmask_m(fix(i),fjx(i))*Wu(fix(i),fjx(i),t)*30*86400/100;
    end
    for i = 1:length(fiy)
        Whi_int = nanmean([Whi(fiy(i)+4,fjy(i)+1,t);Whi(fiy(i)+4,fjy(i)+2,t);Whi(fiy(i)+5,fjy(i)+1,t);Whi(fiy(i)+5,fjy(i)+2,t)]);
        WSIC_int = nanmean([WSIC(fiy(i)+4,fjy(i)+1,t);WSIC(fiy(i)+4,fjy(i)+2,t);WSIC(fiy(i)+5,fjy(i)+1,t);WSIC(fiy(i)+5,fjy(i)+2,t)]);
        volumeTy(i) = Whi_int*ymask_m(fiy(i),fjy(i))*Wv(fiy(i),fjy(i),t)*30*86400/100;
        areaTy(i) = WSIC_int*ymask_m(fiy(i),fjy(i))*Wv(fiy(i),fjy(i),t)*30*86400/100;
    end
    Vimport_monthly_RCP85(n,t,k) = nansum(volumeTx) + nansum(volumeTy);
    Aimport_monthly_RCP85(n,t,k) = nansum(areaTx) + nansum(areaTy);
	clearvars volumeTx volumeTy areaTx areaTy
    end
end
fprintf('finished RCP8.5 ensemble # %d\n',k)
end
Vimport_20C = squeeze(nansum(Vimport_monthly_20C,2));
Aimport_20C = squeeze(nansum(Aimport_monthly_20C,2));
Vimport_RCP85 = squeeze(nansum(Vimport_monthly_RCP85,2));
Aimport_RCP85 = squeeze(nansum(Aimport_monthly_RCP85,2));

save /home/ocean_personal_data/samc/CESM/variables/CESM_ice_import.mat Vimport_monthly_20C Vimport_monthly_RCP85 Aimport_monthly_20C Aimport_monthly_RCP85 Vimport_20C Vimport_RCP85 Aimport_20C Aimport_RCP85

