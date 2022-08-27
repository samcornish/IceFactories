% comparing observational data to CESM -- all 40 models
% ERA5 2m T
% SIC from NASA Goddard
% div from Polar Pathfinder

clear
close all

addpath ~/Documents/MATLAB/project/functions/
addpath ~/Documents/MATLAB/project/functions/export_fig/
addpath ~/Documents/MATLAB/project/functions/cmocean/
addpath ~/Documents/MATLAB/project/data/CESM/
addpath ~/Documents/MATLAB/project/variables/
addpath ~/Documents/MATLAB/project/data/ERA5/
addpath ~/Documents/MATLAB/project/data/sea_ice_index/
addpath ~/Documents/MATLAB/project/data/NOAA_OI_SSTv2/
addpath ~/Documents/MATLAB/project/data/SnowLG/
addpath ~/Documents/MATLAB/project/data/HADISST/

load KL_IceFactories_2mT.mat
X_2mT = X; lat_2mT = lat; lon_2mT = lon; clearvars X lat lon
load KL_Ice_Factories_SIC.mat
X_SIC = X; lat_SIC = lat; lon_SIC = lon; clearvars X lat lon
load PP_div_KL_ts_map.mat
X_PP = X; lat_PP = lat; lon_PP = lon; clearvars X lat lon
load KL_Ice_Factories_SST.mat
X_NOAA = X; lat_SST = lat; lon_SST = lon; Sep_SST_NOAA = Sep_SST_KL; Sep_SST_map_NOAA = Sep_SST_map; clearvars X lat lon Sep_SST_KL Sep_SST_map
load snowLG_KL.mat
load KL_Ice_Factories_HadISST.mat X Sep_SST_KL Sep_SST_map lat lon
X_HadISST = X; lat_HadISST = lat; lon_HadISST = lon; Sep_SST_HadISST = Sep_SST_KL; Sep_SST_map_HadISST = Sep_SST_map; clearvars X lat lon Sep_SST_KL Sep_SST_map


load CESM_T10m_20C_RCP85_all40.mat monthly_T10m_20C monthly_T10m_RCP85
load CESM_iceprod_vars_40.mat

total_area = 1.2543e12;
Sep_open_water_20C = total_area - SepSIA_20C;
Sep_open_water_RCP85 = total_area - SepSIA_RCP85;

deltaT_20C = -(SAT_20C+1.8);
deltaT_RCP85 = -(SAT_RCP85+1.8);

Sep_10mT_20C = squeeze(monthly_T10m_20C(:,1,:));
Sep_10mT_RCP85 = squeeze(monthly_T10m_RCP85(:,1,:));

yrs = 1921:2006;
yrs85 = 2007:2079;
yrs_full = 1921:2079;

comp_areadiv_20C = pos_areadiv_20C - net_areadiv_20C;
comp_areadiv_RCP85 = pos_areadiv_RCP85 - net_areadiv_RCP85;

compdiv_area = posdiv_area - netdiv_area;

% KL outline
vec1lat = linspace(82,82,41);
vec1lon = linspace(70,110,41);
vec2lat = linspace(78,78,41);
vec2lon = linspace(110,150,41);

shape_lat = [72.75,vec1lat,82,78,vec2lat,78,72];
shape_lon = [70,vec1lon,110,110,vec2lon,150,150];

% some colours
c1 = [0 119 187]/256;
c2 = [51 187 238]/256;
c3 = [0 153 136]/256;
c4 = [238 119 51]/256;
c5 = [204 51 17]/256;
c6 = [238 51 119]/256;
c7 = [187 187 187]/256;

c1t = [0 119 187 128]/256;
c2t = [51 187 238 128]/256;
c3t = [0 153 136 128]/256;
c4t = [238 119 51 128]/256;
c5t = [204 51 17 128]/256;
c6t = [238 51 119 128]/256;
c7t = [187 187 187 128]/256;

c_blue = [0.6 0.8 0.95];
c_orange = [0.95 0.8 0.6];

c_bluet = [0.6 0.8 0.95 0.5];
c_oranget = [0.95 0.8 0.6 0.5];
%% create full timeseries

SAT_full = cat(1,SAT_20C,SAT_RCP85);
mSAT_full = mean(SAT_full,2);

deltaT_full = cat(1,deltaT_20C,deltaT_RCP85);
mdeltaT_full = mean(deltaT_full,2);

SepSIA_full = cat(1,SepSIA_20C,SepSIA_RCP85);
mSepSIA_full = mean(SepSIA_full,2);

Sep_open_water_full = cat(1,Sep_open_water_20C,Sep_open_water_RCP85);
mSep_open_water_full = mean(Sep_open_water_full,2);

total_areadiv_full = cat(1,pos_areadiv_20C,pos_areadiv_RCP85);
mtotal_areadiv_full = mean(total_areadiv_full,2);

net_areadiv_full = cat(1,net_areadiv_20C,net_areadiv_RCP85);
mnet_areadiv_full = mean(net_areadiv_full,2);

comp_areadiv_full = total_areadiv_full - net_areadiv_full;
mcomp_areadiv_full = mean(comp_areadiv_full,2);

Sep_10mT_full = cat(1,Sep_10mT_20C,Sep_10mT_RCP85);
mSep_10mT_full = mean(Sep_10mT_full,2);

hs_full = cat(1,hs_20C,hs_RCP85);
mhs_full = mean(hs_full,2);


%% 2m T plot

figure;
ax1 = axes('Position',[0.125,0.11,0.8,0.8]);
p = plot(X_2mT,winter_2mT-273.15); p.LineWidth = 2;
p.Color = c5; hold on;

p = plot(yrs,mean(SAT_20C,2)); p.Color = c_blue; p.LineWidth = 2;
p = plot(yrs85,mean(SAT_RCP85,2)); p.Color = c_orange; p.LineWidth = 2;
L = legend('ERA5','CESM 20C','CESM RCP8.5'); L.Location = 'NorthEast'; L.AutoUpdate = 'off';
p = plot(yrs,SAT_20C,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs,mean(SAT_20C,2)); p.Color = c_blue/2; p.LineWidth = 2;
p = plot(yrs85,SAT_RCP85,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85,mean(SAT_RCP85,2)); p.Color = c_orange/2; p.LineWidth = 2;
p = plot(X_2mT,winter_2mT-273.15); p.LineWidth = 2;
p.Color = c5;
grid on;
ylabel('2 m temperature, ^oC'); %xlim([1978 2020]); 
ylim([-23 -14]);
title('Winter 2 m temperature from ERA5');
ax1.FontName = 'Avenir Next'; ax1.FontSize = 13;
%ax1.XLim = [1978 2020]; 
ax1.YLim = [-28 -10]; grid on;  

ax2 = axes('Position',[0.12,0.45,0.4,0.35]);
hw = worldmap([67.5 90],[0 180]); % plot axes, lat 60-90N
hold on;
hw.View = [270,90];
hw.GridLineStyle = 'none';
set(findall(hw,'Tag','PLabel'),'visible','off')
set(findall(hw,'Tag','MLabel'),'visible','off')
pcolorm(lat_2mT,lon_2mT,(nanmean(T_map(:,:,end-20:end-1),3) - nanmean(T_map(:,:,1:20),3))')
load coastlines
geoshow(coastlat, coastlon, 'Color', 'black');
% geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
g = geoshow(shape_lat,shape_lon,'Color',c_blue); g.LineWidth = 1.5; 
colormap(cmocean('balance','pivot',0)); 
c = colorbar('location','southoutside');
c.AxisLocation = 'out'; c.Position = c.Position + [0.1,-0.075,-0.2,0.02];
c.Label.String = 'change in 2 m temp.'; 
title('1980-1999 to 2000-2019');
ax2.FontName = 'Avenir Next'; ax2.FontSize = 13;

% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/compare_2mT_ERA5_CESM.fig';
% saveas(gcf,fn);

%% SIC

figure;
ax1 = axes('Position',[0.125,0.11,0.8,0.8]);
X = 1979:2019;
p = plot(X,SIA_sum/1e9); p.LineWidth = 2; hold on; p.Color = c5;
p = plot(yrs-1,mean(SepSIA_20C,2)/1e9,'Color',c_blue,'LineWidth',2);
p = plot(yrs85-1,mean(SepSIA_RCP85,2)/1e9,'Color',c_orange,'LineWidth',2);
L = legend('ERA5','CESM 20C','CESM RCP8.5'); L.Location = 'SouthWest'; L.AutoUpdate = 'off';
p = plot(yrs-1,SepSIA_20C/1e9,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs-1,mean(SepSIA_20C,2)/1e9,'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85-1,SepSIA_RCP85/1e9,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85-1,mean(SepSIA_RCP85,2)/1e9,'Color',c_orange/2,'LineWidth',2);
p = plot(X,SIA_sum/1e9); p.LineWidth = 2; hold on; p.Color = c5;

grid on;
ylabel('Sea ice area, 1000 km^2'); ylim([0 1100]); %xlim([1978 2020]);
title('September sea ice area in Kara and Laptev seas');
ax1.FontName = 'Avenir Next'; ax1.FontSize = 13;


ax2 = axes('Position',[0.59,0.475,0.36,0.36]);
hw = worldmap([67.5 90],[0 180]); % plot axes, lat 60-90N
hold on;
hw.View = [270,90];
hw.GridLineStyle = 'none';
set(findall(hw,'Tag','PLabel'),'visible','off')
set(findall(hw,'Tag','MLabel'),'visible','off')
pcolorm(lat_SIC,lon_SIC,mean(SIC_map(:,:,end-19:end),3) - mean(SIC_map(:,:,2:21),3))
load coastlines
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
g = geoshow(shape_lat,shape_lon,'Color',c_orange); g.LineWidth = 1.5; 
colormap(cmocean('balance','pivot',0)); 
c = colorbar('location','southoutside');
c.AxisLocation = 'out'; c.Position = c.Position + [0.056,-0.1,-0.12,0.02];
c.Label.String = 'change in sea ice conc.'; 
title('1980-1999 to 2000-2019');
ax2.FontName = 'Avenir Next'; ax2.FontSize = 13;
c.Label.BackgroundColor = [1 1 1];

% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/compare_SIC_NSIDC_CESM.fig';
% saveas(gcf,fn);

%% div

figure;
ax1 = axes('Position',[0.125,0.1,0.8,0.8]);
p = plot(yrs,mean(pos_areadiv_20C,2),'Color',c2,'LineWidth',2); hold on;
p = plot(yrs,mean(net_areadiv_20C,2),'Color',c7,'LineWidth',2);
p = plot(X_PP,posdiv_area); p.LineWidth = 2; p.Color = c4; hold on;
p = plot(X_PP,netdiv_area); p.LineWidth = 2; p.Color = c6;
grid on;
ylabel('winter divergence, m^2'); L = legend('CESM pos div','CESM net div','PP pos div','PP pos div'); L.Location = 'NorthEast'; L.AutoUpdate = 'off';
p = plot(yrs,pos_areadiv_20C,'Color',c2t,'LineWidth',0.5); hold on;
p = plot(yrs,net_areadiv_20C,'Color',c7t,'LineWidth',0.5);
p = plot(yrs85,pos_areadiv_RCP85,'Color',c2t,'LineWidth',0.5);
p = plot(yrs85,net_areadiv_RCP85,'Color',c7t,'LineWidth',0.5);
p = plot(yrs,mean(pos_areadiv_20C,2),'Color',c2/2,'LineWidth',2); hold on;
p = plot(yrs,mean(net_areadiv_20C,2),'Color',c7/2,'LineWidth',2);
p = plot(yrs85,mean(pos_areadiv_RCP85,2),'Color',c2/2,'LineWidth',2);
p = plot(yrs85,mean(net_areadiv_RCP85,2),'Color',c7/2,'LineWidth',2);
p = plot(X_PP,posdiv_area); p.LineWidth = 2; p.Color = c4; hold on;
p = plot(X_PP,netdiv_area); p.LineWidth = 2; p.Color = c6;

t1 = text(1990, 0.85,'20C'); t1.FontName = 'Avenir Next'; t1.FontSize = 13;
t2 = text(2010, 0.85,'RCP8.5'); t2.FontName = 'Avenir Next'; t2.FontSize = 13;

title('Winter divergence from Polar Pathfinder');
ax1.FontName = 'Avenir Next'; ax1.FontSize = 13; %ax1.YLim = [0 1.2]; %ax1.XLim = [1975 2020];

ax2 = axes('Position',[0.135,0.575,0.36,0.36]);
hw = worldmap([70 85],[60 160]); % plot axes
hold on;
% hw.View = [270,90];
hw.GridLineStyle = 'none';
set(findall(hw,'Tag','PLabel'),'visible','off')
set(findall(hw,'Tag','MLabel'),'visible','off')
pcolorm(lat_PP,lon_PP,nanmean(Wdiv_map(:,:,end-18:end),3) - nanmean(Wdiv_map(:,:,1:20),3));
%pcolorm(lat_PP,lon_PP,nanmean(Wdiv_map,3))
load coastlines
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
g = geoshow(shape_lat,shape_lon,'Color',c3); g.LineWidth = 1.5; 
colormap(cmocean('balance','pivot',0)); 
c = colorbar('location','southoutside');
c.AxisLocation = 'out'; c.Position = c.Position + [0.056,-0.05,-0.12,0.02];
c.Label.String = 'change in div, %/day'; c.Label.BackgroundColor = [1 1 1];
% title('First vs last decade');
ax2.FontName = 'Avenir Next'; ax2.FontSize = 13;
t3= text(ax1,1980,2.6e12,'1980-1999 to \newline     2000-2018'); t3.FontName = 'Avenir Next'; t3.FontSize = 13; t3.FontWeight = 'bold'; 

% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/compare_div_PP_CESM.fig';
% saveas(gcf,fn);

%% September SST

figure;
ax1 = axes('Position',[0.125,0.11,0.8,0.8]);
p = plot(X_NOAA,Sep_SST_NOAA); p.LineWidth = 2; hold on; p.Color = c5;
p = plot(X_HadISST,Sep_SST_HadISST); p.LineWidth = 2; hold on; p.Color = c4;
p = plot(yrs-1,mean(Sep_10mT_20C,2),'Color',c_blue,'LineWidth',2);
p = plot(yrs85-1,mean(Sep_10mT_RCP85,2),'Color',c_orange,'LineWidth',2);
L = legend('SST - NOAA OI SSTv2','SST - HadISST','10 m T - CESM 20C','10 m T - CESM RCP8.5'); L.Location = 'SouthEast'; L.AutoUpdate = 'off';
p = plot(yrs-1,Sep_10mT_20C,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs-1,mean(Sep_10mT_20C,2),'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85-1,Sep_10mT_RCP85,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85-1,mean(Sep_10mT_RCP85,2),'Color',c_orange/2,'LineWidth',2);
p = plot(X_NOAA,Sep_SST_NOAA); p.LineWidth = 2; hold on; p.Color = c5;

grid on;
ylabel('^oC'); ylim([-2 6]); %xlim([1978 2022]);
title('September SST in Kara and Laptev seas from NOAA OI SSTv2');
ax1.FontName = 'Avenir Next'; ax1.FontSize = 13;

ax2 = axes('Position',[0.12,0.48,0.4,0.35]);
hw = worldmap([67.5 90],[0 180]); % plot axes, lat 60-90N
hold on;
hw.View = [270,90];
hw.GridLineStyle = 'none';
set(findall(hw,'Tag','PLabel'),'visible','off')
set(findall(hw,'Tag','MLabel'),'visible','off')
pcolorm(lat_HadISST,lon_HadISST,(nanmean(Sep_SST_map_HadISST(:,:,40:end),3) - nanmean(Sep_SST_map_HadISST(:,:,20:39),3)))
load coastlines
%geoshow(coastlat, coastlon, 'Color', 'black');
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
g = geoshow(shape_lat,shape_lon,'Color',c2); g.LineWidth = 1.5; 
colormap(cmocean('balance','pivot',0)); 
c = colorbar('location','southoutside');
c.AxisLocation = 'out'; c.Position = c.Position + [0.1,-0.1,-0.2,0.02];
c.Label.String = 'change in SST, ^oC'; 
title('1980-1999 to 2000-2020');
ax2.FontName = 'Avenir Next'; ax2.FontSize = 13;

% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/compare_SST_NOAA_HadISST_CESM.fig';
% saveas(gcf,fn);


%% Just the timeseries

M = subplot2axes(6,1,0.1,0.05); %M(:,1) = M(:,1) + 0.035;
figure;
% SAT
ax1 = axes('Position',M(1,:));
p = plot(X_2mT,winter_2mT-273.15); p.LineWidth = 2;
L = legend('ERA5'); L.Location = 'West'; L.AutoUpdate = 'off';
p.Color = c5; hold on;
p = plot(yrs,SAT_20C,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs,mean(SAT_20C,2)); p.Color = c_blue/2; p.LineWidth = 2;
p = plot(yrs85,SAT_RCP85,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85,mean(SAT_RCP85,2)); p.Color = c_orange/2; p.LineWidth = 2;
p = plot(X_2mT,winter_2mT-273.15); p.LineWidth = 2;
p.Color = c5;
grid on;
ylabel('^oC'); %xlim([1978 2020]); 
ax1.YLim = [-30 -10]; grid on;  
yyaxis right
p1 = plot(yrs_full,mdeltaT_full); ylabel('\DeltaT, ^oC'); p1.LineWidth = 1;
title('Mean winter SAT');
ax1.FontName = 'Helvetica'; ax1.FontSize = 13;
%ax1.XLim = [1978 2020]; 
ax1.YLim = [10 30]; grid on;  

%SIC
ax2 = axes('Position',M(2,:));
p = plot(X_SIC,SIA_sum/1e12); p.LineWidth = 2; hold on; p.Color = c5;
p1 = plot(yrs_full-1,mSep_open_water_full/1e12); p1.LineWidth = 1; p1.Color = [0 0.4470 0.7410];
L = legend('NASA Goddard','Sep open water'); L.Location = 'East'; L.AutoUpdate = 'off';
p = plot(yrs-1,SepSIA_20C/1e12,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs-1,mean(SepSIA_20C,2)/1e12,'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85-1,SepSIA_RCP85/1e12,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85-1,mean(SepSIA_RCP85,2)/1e12,'Color',c_orange/2,'LineWidth',2);
p = plot(X_SIC,SIA_sum/1e12); p.LineWidth = 2; hold on; p.Color = c5;
p1 = plot(yrs_full-1,mSep_open_water_full/1e12); p1.LineWidth = 1; p1.Color = [0 0.4470 0.7410];
grid on;
ylabel('10^6 km^2'); ylim([0 1.26]); %xlim([1978 2020]);
title('Sep SIA');
%ax1.XLim = [1978 2020]; 
ylim([0 1.26]); grid on; 
ax2.FontName = 'Helvetica'; ax2.FontSize = 13;

% div
ax3 = axes('Position',M(3,:));
p = plot(X_PP,netdiv_area/1e12); p.LineWidth = 2; p.Color = c5; hold on;
L = legend('Polar Pathfinder'); L.Location = 'NorthWest'; L.AutoUpdate = 'off';
p = plot(yrs,mean(net_areadiv_20C,2)/1e12,'Color',c_blue,'LineWidth',2);
grid on;
ylabel('10^6 km^2'); %L = legend('CESM pos div','CESM net div','PP pos div','PP pos div'); L.Location = 'best'; L.AutoUpdate = 'off';
p = plot(yrs,net_areadiv_20C/1e12,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs85,net_areadiv_RCP85/1e12,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs,mean(net_areadiv_20C,2)/1e12,'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85,mean(net_areadiv_RCP85,2)/1e12,'Color',c_orange/2,'LineWidth',2);
p = plot(X_PP,netdiv_area/1e12); p.LineWidth = 2; p.Color = c5;
title('Net winter sea ice area diverged');
ax3.FontName = 'Helvetica'; ax3.FontSize = 13; ax3.YLim = [0 2]; %ax1.XLim = [1975 2020];


% comp div 
ax4 = axes('Position',M(4,:));
p = plot(X_PP,compdiv_area/1e12); p.LineWidth = 2; p.Color = c5; hold on;
L = legend('Polar Pathfinder'); L.Location = 'NorthWest'; L.AutoUpdate = 'off';
p = plot(yrs,mean(comp_areadiv_20C,2)/1e12,'Color',c_blue,'LineWidth',2); 
grid on;
ylabel('10^6 km^2'); %L = legend('CESM pos div','CESM net div','PP pos div','PP pos div'); L.Location = 'best'; L.AutoUpdate = 'off';
p = plot(yrs,comp_areadiv_20C/1e12,'Color',c_bluet,'LineWidth',0.5); 
p = plot(yrs85,comp_areadiv_RCP85/1e12,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs,mean(comp_areadiv_20C,2)/1e12,'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85,mean(comp_areadiv_RCP85,2)/1e12,'Color',c_orange/2,'LineWidth',2);
p = plot(X_PP,compdiv_area/1e12); p.LineWidth = 2; p.Color = c5; hold on;
title('Compensated winter sea ice area diverged');
ax4.FontName = 'Helvetica'; ax4.FontSize = 13; ax4.YLim = [0 2]; %ax1.XLim = [1975 2020];


% SST
ax5 = axes('Position',M(5,:));
p = plot(X_HadISST,Sep_SST_HadISST); p.LineWidth = 2; hold on; p.Color = c4;
p = plot(X_NOAA,Sep_SST_NOAA); p.LineWidth = 2; hold on; p.Color = c6;
L = legend('NOAA','HadISST'); L.Location = 'NorthWest'; L.AutoUpdate = 'off';
p = plot(yrs-1,mean(Sep_10mT_20C,2),'Color',c_blue,'LineWidth',2);
p = plot(yrs85-1,mean(Sep_10mT_RCP85,2),'Color',c_orange,'LineWidth',2);
% L = legend('SST - NOAA OI SSTv2','10 m T - CESM 20C','10 m T - CESM RCP8.5'); L.Location = 'SouthEast'; L.AutoUpdate = 'off';
p = plot(yrs-1,Sep_10mT_20C,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs-1,mean(Sep_10mT_20C,2),'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85-1,Sep_10mT_RCP85,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85-1,mean(Sep_10mT_RCP85,2),'Color',c_orange/2,'LineWidth',2);
p = plot(X_HadISST,Sep_SST_HadISST); p.LineWidth = 2; hold on; p.Color = c4;
p = plot(X_NOAA,Sep_SST_NOAA); p.LineWidth = 2; hold on; p.Color = c6;
grid on;
ylabel('^oC'); ylim([-2 6]); %xlim([1978 2022]);
title('Sep SST');
ax5.FontName = 'Helvetica'; ax5.FontSize = 13;

% snow
ax6 = axes('Position',M(6,:));
p = plot(snowKLM2.yrs,snowKLM2.Mean*100); p.LineWidth = 2; hold on; p.Color = c5;
L = legend('SnowModel-LG'); L.Location = 'NorthWest'; L.AutoUpdate = 'off';
p = plot(yrs,mean(hs_20C,2)*100,'Color',c_blue,'LineWidth',2);
p = plot(yrs85,mean(hs_RCP85,2)*100,'Color',c_orange,'LineWidth',2);
% L = legend('SST - NOAA OI SSTv2','10 m T - CESM 20C','10 m T - CESM RCP8.5'); L.Location = 'SouthEast'; L.AutoUpdate = 'off';
p = plot(yrs,hs_20C*100,'Color',c_bluet,'LineWidth',0.5);
p = plot(yrs,mean(hs_20C,2)*100,'Color',c_blue/2,'LineWidth',2);
p = plot(yrs85,hs_RCP85*100,'Color',c_oranget,'LineWidth',0.5);
p = plot(yrs85,mean(hs_RCP85,2)*100,'Color',c_orange/2,'LineWidth',2);
p = plot(snowKLM2.yrs,snowKLM2.Mean*100); p.LineWidth = 2; hold on; p.Color = c5;
grid on;
ylabel('cm'); ylim([0 30]); %xlim([1978 2022]);
title('Mean winter snow depth');
ax6.FontName = 'Helvetica'; ax6.FontSize = 13;
ax1.XTick = [1940 1960 1980 2000 2020 2040 2060]; ax2.XTick = [1940 1960 1980 2000 2020 2040 2060]; ax3.XTick = [1940 1960 1980 2000 2020 2040 2060]; 
ax4.XTick = [1940 1960 1980 2000 2020 2040 2060]; ax5.XTick = [1940 1960 1980 2000 2020 2040 2060]; ax6.XTick = [1940 1960 1980 2000 2020 2040 2060];
ax1.XTickLabel = []; ax2.XTickLabel = []; ax3.XTickLabel = []; ax4.XTickLabel = []; ax5.XTickLabel = []; %ax6.XTickLabel = {'1940','','1980','', '2020','', '2060'};
ax1.XLim = [1920 2080]; ax2.XLim = [1920 2080]; ax3.XLim = [1920 2080]; ax4.XLim = [1920 2080]; ax5.XLim = [1920 2080]; ax6.XLim = [1920 2080];


%%
fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/compare_timeseries_CESM_obsv6_40.fig';
saveas(gcf,fn);
% % % 
set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/compare_timeseries_CESM_obsv6_40.png -m3

print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/compare_timeseries_CESM_obsv7_40.pdf','-dpdf','-bestfit');
