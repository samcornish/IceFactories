% mask of LAPTEV and KARA seas for CESM

clear
close all

addpath /home/ocean_shared_data2/CESM_arctic/divu/

load coastlines

lat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','TLAT');
lon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','TLON');

mask = ones(291,48);    % mask of ones with same size as lat,lon

fE = find(lon>150);
mask(fE) = NaN;
fW = find(lon<70);
mask(fW) = NaN;
fN = find(lat>82);
mask(fN) = NaN;
fL1 = find(lat>78); fL2 = find(lon>110);
fL = intersect(fL1,fL2);
mask(fL) = NaN;

div = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.divu_nh.192001-200512.nc','divu');

d11 = mask.*div(:,:,11);
% I also need a map for the higher resolution runs
lat1xx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.divu_nh.192001-200512.nc','TLAT');
lon1xx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.divu_nh.192001-200512.nc','TLON');

mask1xx = ones(size(lat1xx)); %320*49
fE = find(lon1xx>150);
mask1xx(fE) = NaN;
fW = find(lon1xx<70);
mask1xx(fW) = NaN;
fN = find(lat1xx>82);
mask1xx(fN) = NaN;
fL1 = find(lat1xx>78); fL2 = find(lon1xx>110);
fL = intersect(fL1,fL2);
mask1xx(fL) = NaN;

dxt1xx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.divu_nh.192001-200512.nc','dxt');
dyt1xx = ncread('b.e11.B20TRC5CNBDRD.f09_g16.101.cice.h.divu_nh.192001-200512.nc','dyt');
cellarea1xx = dxt1xx.*dyt1xx;
%figure;
%hw = worldmap([60 90],[0 180]); % plot axes, lat 60-90N
%hw.View = [270,90];
%hw.GridLineStyle = 'none';
%set(findall(hw,'Tag','PLabel'),'visible','off')
%set(findall(hw,'Tag','MLabel'),'visible','off')
%pcolorm(double(lat),double(lon),d11)
%geoshow(coastlat, coastlon, 'Color','black');

save /home/ocean_personal_data/samc/CESM/variables/KLmask.mat mask lat lon mask1xx lat1xx lon1xx dxt1xx dyt1xx cellarea1xx
%% mask for velocities
addpath /home/ocean_shared_data2/CESM_arctic/uvel/

Ulat = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.uvel_nh.192001-200512.nc','ULAT');
Ulon = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.uvel_nh.192001-200512.nc','ULON');

Umask = ones(size(Ulat));    % mask of ones with same size as lat,lon

fE = find(Ulon>150);
Umask(fE) = NaN;
fW = find(Ulon<70);
Umask(fW) = NaN;
fN = find(Ulat>82);
Umask(fN) = NaN;
fL1 = find(Ulat>78); fL2 = find(Ulon>110);
fL = intersect(fL1,fL2);
Umask(fL) = NaN;

%v = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.vvel_nh.192001-200512.nc','vvel');
%u = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.uvel_nh.192001-200512.nc','uvel');
%v11 = Umask.*mean(v,3);
%u11 = Umask.*mean(u,3);
%figure;
%hw = worldmap([70 85],[60 160]); % plot axes, lat 60-90N
% hw.View = [270,90];
%hw.GridLineStyle = 'none';
%set(findall(hw,'Tag','PLabel'),'visible','off')
%set(findall(hw,'Tag','MLabel'),'visible','off')
%pcolorm(double(Ulat),double(Ulon),u11)
%geoshow(coastlat, coastlon, 'Color','black');

save /home/ocean_personal_data/samc/CESM/variables/KLmask.mat Umask Ulat Ulon -append

%% export masks
u = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.uvel_nh.192001-200512.nc','uvel');
dxu = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.uvel_nh.192001-200512.nc','dxu');
dyu = ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.uvel_nh.192001-200512.nc','dyu');
um = mean(u,3);
Umask(isnan(Umask)) = 0;
[fi fj] = find(Umask);
Umask(isnan(um)) = NaN;
fyn = [1,1]; fys = [1,1];
fxe = [1,1]; fxw = [1,1];
ymask = zeros(size(Umask));
xmask = zeros(size(Umask));

for i = 1:length(fi)
if Umask(fi(i),fj(i)+1)==0
    fyn = cat(1,fyn,[fi(i),fj(i)]);
    ymask(fi(i),fj(i)) = -1;
end
if Umask(fi(i),fj(i)-1)==0
    fys = cat(1,fys,[fi(i),fj(i)]);
    ymask(fi(i),fj(i)) = 1;
end
if Umask(fi(i)+1,fj(i))==0
    fxe = cat(1,fxe,[fi(i),fj(i)]);
    xmask(fi(i),fj(i)) = -1;
end
if Umask(fi(i)-1,fj(i))==0
    fxw = cat(1,fxw,[fi(i),fj(i)]);
    xmask(fi(i),fj(i)) = 1;
end
end
save /home/ocean_personal_data/samc/CESM/variables/KLmask.mat xmask ymask dxu dyu -append
