% script to plot decomposition of ice production into freezing area days
% and growth rate

clear
close all

addpath ~/Documents/MATLAB/project/functions/
addpath ~/Documents/MATLAB/project/functions/addpath_recurse/
addpath ~/Documents/MATLAB/project/functions/export_fig/
addpath ~/Documents/MATLAB/project/functions/cmocean/
addpath ~/Documents/MATLAB/project/data/CESM/
addpath ~/Documents/MATLAB/project/variables/
addpath_recurse ~/Documents/MATLAB/project/functions/boundedline-pkg/

yrs = 1921:2006;
yrs85 = 2007:2079;
yrs_full = 1921:2079;

load CESM_iceprod_vars_40.mat


miceprod_20C = mean(iceprod_20C,2);
miceprod_RCP85 = mean(iceprod_RCP85,2);
mfrz_area_days_20C = mean(frz_area_days_20C,2);
mfrz_area_days_RCP85 = mean(frz_area_days_RCP85,2);
growth_rate_20C = iceprod_20C./frz_area_days_20C;
growth_rate_RCP85 = iceprod_RCP85./frz_area_days_RCP85;
mgrowth_rate_20C = mean(growth_rate_20C,2);
mgrowth_rate_RCP85 = mean(growth_rate_RCP85,2);
%% figure: decomposing ice production into growth rate and freezing area days
yrs = 1921:2006;

c_blue = [0.6 0.8 0.95];
c_orange = [0.95 0.8 0.6];

c_bluet = [0.6 0.8 0.95 0.5];
c_oranget = [0.95 0.8 0.6 0.5];

figure;
ax1 = axes('Position',[0.32,0.6,0.36,0.35]);
p = plot(yrs,-iceprod_20C/1e12,'Color',c_bluet); hold on;
p = plot(yrs,-miceprod_20C/1e12); p.LineWidth = 2; p.Color = c_blue/2;
p = plot(yrs85,-iceprod_RCP85/1e12,'Color',c_oranget); 
p = plot(yrs85,-miceprod_RCP85/1e12); p.LineWidth = 2; p.Color = c_orange/2;
title('Total winter ice production'); ax1.XLim = [1920 2080]; ylabel('x1000 km^3');
t1 = text(1980,1.2,'20C'); t1.Color = c_blue; t2 = text(2020,1.2,'RCP85'); t2.Color = c_orange; 
t1.FontSize = 13; t2.FontSize = 13; t1.FontName = 'Helvetica'; t2.FontName = 'Helvetica';
ta = text(ax1,1905,2.75,'a)'); ta.FontName = 'Helvetica'; ta.FontSize = 13; ta.FontWeight = 'bold';
ax1.FontSize = 14; ax1.FontName = 'Helvetica'; grid on;
ax1.XTick = [1920:20:2080]; ax1.XTickLabel = {'1920','','1960','','2000','','2040','','2080'};

total_area = 1.2543e12;

%ax3 is first because I'm reordering the plots and lazy
ax3 = axes('Position',[0.12,0.1,0.36,0.35]);
p = plot(yrs,frz_area_days_20C/total_area,'Color',c_bluet); hold on;
p = plot(yrs,mfrz_area_days_20C/total_area); p.LineWidth = 2; p.Color = c_blue/2;
p = plot(yrs85,frz_area_days_RCP85/total_area,'Color',c_oranget);
p = plot(yrs85,mfrz_area_days_RCP85/total_area); p.LineWidth = 2; p.Color = c_orange/2;
title('Freezing area days, {\it "Usage"}'); ax3.XLim = [1920 2080]; ylabel('x 1.25x10^6 km^2 day'); ax3.YAxisLocation = 'Left';
tb = text(ax3,1905,230,'b)'); tb.FontName = 'Helvetica'; tb.FontSize = 13; tb.FontWeight = 'bold';
ax3.FontSize = 14; ax3.FontName = 'Helvetica'; grid on; 
ax3.XTick = [1920:20:2080]; ax3.XTickLabel = {'1920','','1960','','2000','','2040','',''};


ax2 = axes('Position',[0.53,0.1,0.36,0.35]);
p = plot(yrs,-growth_rate_20C*1e3,'Color',c_bluet); hold on;
p = plot(yrs,-mgrowth_rate_20C*1e3); p.LineWidth = 2; p.Color = c_blue/2;
p = plot(yrs85,-growth_rate_RCP85*1e3,'Color',c_oranget);
p = plot(yrs85,-mgrowth_rate_RCP85*1e3); p.LineWidth = 2; p.Color = c_orange/2;
title('Growth rate, {\it "Efficiency"}'); ax2.XLim = [1920 2080]; ylabel('mm/day');
tc = text(ax2,1905,11,'c)'); tc.FontName = 'Helvetica'; tc.FontSize = 13; tc.FontWeight = 'bold';
ax2.FontSize = 14; ax2.FontName = 'Helvetica'; grid on;
ax2.XTick = [1920:20:2080]; ax2.XTickLabel = {'1920','','1960','','2000','','2040','','2080'}; ax2.YAxisLocation = 'Right';



% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/IceProd_decomp_power_usage.fig';
% saveas(gcf,fn);

set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/IceProd_decomp_power_usage.pdf 
print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/IceProd_decomp_power_usage2_40.pdf','-dpdf','-bestfit')