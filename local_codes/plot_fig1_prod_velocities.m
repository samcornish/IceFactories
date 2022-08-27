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


load CESM_spatial_maps.mat
load bathymetry_Arctic_ocean.mat bathymetry

lat = double(ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','TLAT'));
lon = double(ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.cice.h.hi_nh.192001-200512.nc','TLON'));

latpop = double(ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.IFRAC.192001-200512.nc','TLAT'));
lonpop = double(ncread('b.e11.B20TRC5CNBDRD.f09_g16.002.pop.h.IFRAC.192001-200512.nc','TLONG'));

yrs = 1921:2006;
yrs85 = 2007:2079;
yrs_full = 1921:2079;

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

%1980-1999
i8099 = 60:79; % apply this to yrs
%2000-2019
% will need to concatenate end of 20C to RCP85
i0019a = 80:86;
i0019b = 1:13;

i7920a = 59:86;
i7920b = 1:14;

%% figure 1a

% 1979-2020 for parity with Polar Pathfinder
M = subplot2axes(1,2,0.05,0.05);

figure;
ax1 = axes('Position',M(1,:));
hw = worldmap([67.5 90],[0 180]); % plot axes, lat 60-90N
hold on;
hw.View = [270,90];
hw.GridLineStyle = 'none';
set(findall(hw,'Tag','PLabel'),'visible','off')
set(findall(hw,'Tag','MLabel'),'visible','off')
pcolorm(lat,lon,-1000*mean(cat(3,I2O_frz_map_20C(:,:,i7920a),I2O_frz_map_RCP85(:,:,i7920b)),3)) % convert to mm/day
% contourm(lat,lon,-1000*mean(I2O_frz_map_20C(:,:,i8099),3),[0:2:20],'k')
load coastlines
%geoshow(coastlat, coastlon, 'Color', 'black');
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
g = geoshow(shape_lat,shape_lon,'Color',c2); g.LineWidth = 1.5; 
% g2 = geoshow(lat1,lon1,'Color',c6); g2.LineWidth = 1.5; g2.LineStyle = '-';
caxis([0 15]);
% colormap(ax1,'hot');
colormap(ax1,cmocean('amp')); 
% contourm(lat,lon,-1000*mean(cat(3,I2O_frz_map_20C(:,:,i7920a),I2O_frz_map_RCP85(:,:,i7920b)),3),[2.5 5 7.5 10 12.5 15],'color','k');
cb1 = colorbar('location','southoutside'); 
cb1.AxisLocation = 'out'; %cb1.Position = cb1.Position + [0.025,-0.065,0,-0.02];
cb1.Position = cb1.Position + [0.1,+0.05,-0.2,0.025];
cb1.Label.String = 'mm/day'; 
% cb1.Ticks = [0:2.5:15]; cb1.TickLabels = {'0','','5','','10','','15'}; cb1.TickLength = 0.2275; 
contourm(latpop,lonpop,bathymetry,[100 500 1000],'color',[1 1 1]);
% cb1.Ticks = [0:2:16]; cb1.TickLength = 0.08;
%c.TickLabels = {'0','','2','','4','','6'}; 
% title('Sea ice production, 1980-1999'); 
title({'Winter ice production, CESM-LE'});
ax1.FontName = 'Helvetica'; ax1.FontSize = 13;

%% figure 1b
clearvars -except M

addpath ~/Documents/MATLAB/project/variables/
addpath ~/Documents/MATLAB/project/functions/cmocean/
addpath ~/Documents/MATLAB/project/data/PolarPathfinder/
addpath ~/Documents/MATLAB/project/functions/
addpath ~/Documents/MATLAB/project/functions/export_fig/

load EASE_mask_polarhole.mat
load coastlines
[ki,kj] = mask_select_EASE(8:15,new_msk);    % function to index area concerned
% KL seas are 17:18 in this EASE mask

load PP_div_polarhole_Gnld1979-2020.mat
load PP_div_polarhole_ArcBas1979-2020.mat u_corr v_corr
speed00 = speedR;
u00 = u_corr;
v00 = v_corr;

speed00KL = zeros(size(speed00))*NaN;
for i = 1:length(ki)
    speed00KL(ki(i),kj(i)) = speed00(ki(i),kj(i));
end
LAT = lat(ki,kj);LAT = LAT(:,1);
LON = lon(ki,kj); LON = LON(:,1);

vec1lat = linspace(82,82,41);
vec1lon = linspace(70,110,41);
vec2lat = linspace(78,78,41);
vec2lon = linspace(110,150,41);

shape_lat = [72.75,vec1lat,82,78,vec2lat,78,72];
shape_lon = [70,vec1lon,110,110,vec2lon,150,150];

ax2 = axes('Position',M(2,:));

k0=find(speed00==0);
speed00(k0) = NaN;

hw = worldmap([67.5 90],[0 180]); % plot axes, lat 60-90N
hw.View = [270,90];
hw.GridLineStyle = 'none';
set(findall(hw,'Tag','PLabel'),'visible','off')
set(findall(hw,'Tag','MLabel'),'visible','off')
pcolorm(lat,lon,speed00); hold on;
%cf = contourfm(lat,lon,div_Tint,21); % plot mesh onto axes
% geoshow(coastlat, coastlon, 'Color', 'black');
quivermc(lat,lon,u00,v00); % ,'color','b'
caxis([0,7])   % for speed
cmap = cmocean('speed');
% caxis([-2e-3, 2e-3])  % for div
% cmap = cmocean('balance','pivot',0);
colormap(ax2,cmap)
g = geoshow(shape_lat,shape_lon,'Color',[0.4 0.6 0.95]); g.LineWidth = 2; 
geoshow('landareas.shp','FaceColor',[0.8 0.8 0.8])
cb2 = colorbar('location','southoutside'); 
cb2.AxisLocation = 'out'; %cb1.Position = cb1.Position + [0.025,-0.065,0,-0.02];
cb2.Position = cb2.Position + [0.1,+0.05,-0.2,0.025];
cb2.Label.String = 'cm/s'; 
title('Mean sea ice velocities, satellite')
ax2.FontSize = 13;
ax2.FontName = 'Helvetica';

% saveas(gcf,'~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/Figure1_a-prod_b-velocitiesv3.fig')

%  axes( 'Position', [0, 0.97, 1, 0.05] ) ;
%  set( gca, 'Color', 'None', 'XColor', 'none', 'YColor', 'none' ) ;
%  text( 0.5, 0, 'Kara and Laptev seas: ice factories of the Arctic Ocean (1979-2020)', 'FontSize', 16', 'FontWeight', 'Bold', ...
%       'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom' ,'FontName','Helvetica') ;
set(gcf, 'Color', 'w');
export_fig '~/Documents/MATLAB/project/data/CESM/figures/Figure1_a-prod_b-velocitiesv4.pdf'
% orient(gcf,'landscape');
% print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/Figure1_a-prod_b-velocitiesv4i.pdf','-dpdf','-bestfit')

