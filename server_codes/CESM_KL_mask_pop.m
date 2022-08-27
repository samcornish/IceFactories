clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/IFRAC/
load coastlines

lat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.IFRAC.192001-200512.nc','TLAT');
lon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.IFRAC.192001-200512.nc','TLONG');

mask = ones(size(lat));    % mask of ones with same size as lat,lon

fE = find(lon>150);
mask(fE) = NaN;
fW = find(lon<70);
mask(fW) = NaN;
fN = find(lat>82);
mask(fN) = NaN;
fL1 = find(lat>78); fL2 = find(lon>110);
fL = intersect(fL1,fL2);
mask(fL) = NaN;

SIC = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.IFRAC.192001-200512.nc','IFRAC');

d11 = mask.*SIC(:,:,11);

%figure;
%hw = worldmap([60 90],[0 180]); % plot axes, lat 60-90N
%hw.View = [270,90];
%hw.GridLineStyle = 'none';
%set(findall(hw,'Tag','PLabel'),'visible','off')
%set(findall(hw,'Tag','MLabel'),'visible','off')
%pcolorm(double(lat),double(lon),d11)
%geoshow(coastlat, coastlon, 'Color','black');

save /home/ocean_personal_data/samc/CESM/variables/KLmask_pop.mat mask lat lon
