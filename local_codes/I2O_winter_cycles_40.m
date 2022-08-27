% plot changing I2O winter cycles--all 40 ens members

clear
close all

addpath ~/Documents/MATLAB/project/functions/
addpath ~/Documents/MATLAB/project/functions/cmocean/
addpath ~/Documents/MATLAB/project/data/CESM/
addpath ~/Documents/MATLAB/project/variables/
addpath ~/Documents/MATLAB/project/functions/export_fig/

load CESM_monthly_I2O_40.mat

yrs = 1921:2006;
yrs85 = 2007:2079;

mmonthly_I2O_20C = mean(monthly_I2O_20C,3);
mmonthly_I2O_RCP85 = mean(monthly_I2O_RCP85,3);
mmonthly_iceprod_20C = mean(monthly_iceprod_20C,3);
mmonthly_iceprod_RCP85 = mean(monthly_iceprod_RCP85,3);
mmonthly_melt_20C = mean(monthly_melt_20C,3);
mmonthly_melt_RCP85 = mean(monthly_melt_RCP85,3);


cmap1 = cmocean('thermal',157);
cmap1t = cat(2,cmap1,0.4*ones(157,1)); % making transparent

% some colours
c1 = [0 119 187]/256;
c2 = [51 187 238]/256;
c3 = [0 153 136]/256;
c4 = [238 119 51]/256;
c5 = [204 51 17]/256;
c6 = [238 51 119]/256;
c7 = [187 187 187]/256;
%% peak timing

% timing of the peak
peak_imb_20C = zeros(84,40,2);
peak_imb_RCP85 = zeros(73,40,2);

for j = 1:40
    for i = 1:86
     [val ind] = max(-monthly_I2O_20C(i,:,j));
     peak_imb_20C(i,j,:) = [val ind];
     [val ind] = max(-monthly_iceprod_20C(i,:,j));
     peak_iceprod_20C(i,j,:) = [val ind];
    end
    for i = 1:73
     [val ind] = max(-monthly_I2O_RCP85(i,:,j));
     peak_imb_RCP85(i,j,:) = [val ind];
     [val ind] = max(-monthly_iceprod_RCP85(i,:,j));
     peak_iceprod_RCP85(i,j,:) = [val ind];
    end
end
mpeak_imb_20C = mean(peak_imb_20C(:,:,2),2);
mpeak_imb_RCP85 = mean(peak_imb_RCP85(:,:,2),2);
mpeak_iceprod_20C = mean(peak_iceprod_20C(:,:,2),2);
mpeak_iceprod_RCP85 = mean(peak_iceprod_RCP85(:,:,2),2);
figure; 
ax1 = axes('Position',[0.125,0.1,0.8,0.8]);
yyaxis left
p1 = plot(yrs,mpeak_imb_20C); p1.LineWidth = 2; hold on; 
p2 = plot(yrs85,mpeak_imb_RCP85); p2.LineWidth = 2; p2.LineStyle = '-';
p = plot(yrs,mpeak_iceprod_20C); p.LineWidth = 2; p.LineStyle = '--';
p = plot(yrs85,mpeak_iceprod_RCP85); p.LineWidth = 2; p.LineStyle = '--';

ax1.YTick = 2:5; 
ax1.YTickLabel = {'Nov','Dec','Jan','Feb'};

yyaxis right
mpeakval_imb_20C = mean(peak_imb_20C(:,:,1),2); 
mpeakval_imb_RCP85 = mean(peak_imb_RCP85(:,:,1),2); 
mpeak_iceprodval_20C = mean(peak_iceprod_20C(:,:,1),2);
mpeak_iceprodval_RCP85 = mean(peak_iceprod_RCP85(:,:,1),2);
p3 = plot(yrs,mpeakval_imb_20C); p3.LineWidth = 2;
p4 = plot(yrs85,mpeakval_imb_RCP85); p4.LineWidth = 2; p4.LineStyle = '-';
p = plot(yrs,mpeak_iceprodval_20C); p.LineWidth = 2; p.LineStyle = '--';
p = plot(yrs85,mpeak_iceprodval_RCP85); p.LineWidth = 2; p.LineStyle = '--';

%% decadal means
decyrs = zeros(16,10);
decm_iceprod = zeros(16,7);
i = 1;
    decm_iceprod(i,:) = mean(mmonthly_iceprod_20C((i-1)*10+1:i*10,:),1);
    decyrs(i,1:9) = yrs((i-1)*10+1:i*10-1);
for i = 2:8
    decm_iceprod(i,:) = mean(mmonthly_iceprod_20C((i-1)*10:i*10-1,:),1);
    decyrs(i,:) = yrs((i-1)*10:i*10-1);
end
    i =9;
    decm_iceprod(i,:) = mean(cat(1,mmonthly_iceprod_20C((i-1)*10+1:(i-1)*10+4,:),mmonthly_iceprod_RCP85(1:3,:)),1);
    decyrs(i,1:7) = cat(2,yrs((i-1)*10+1:(i-1)*10+4),yrs85(1:3));
for i = 10:16
    decm_iceprod(i,:) = mean(mmonthly_iceprod_RCP85((i-10)*10+4:(i-10)*10+13,:),1);
    decyrs(i,:) = yrs85((i-10)*10+4:(i-10)*10+13);
end


cmap2 = cmocean('thermal',16);
cmap2t = [cmap2,0.7*ones(16,1)];
figure;
ax1 = axes('Position',[0.125,0.1,0.8,0.8]);
for i = 1:16
    p = plot(1:7,-decm_iceprod(i,:)/1e9); p.Color = cmap2t(i,:); p.LineWidth = 2;
    hold on;
end
ax1.XTickLabel = {'Oct','Nov','Dec','Jan','Feb','Mar','Apr'};
ax1.YLabel.String = 'Ice production, km^3/month';
% L = legend('1920-1929','1930-1939','1940-1949','1950-1959','1960-1969','1970-1979','1980-1989','1990-1999','2000-2009','2010-2019','2020-2029',...
%     '2030-2039','2040-2049','2050-2059','2060-2069','2070-2079'); L.Location = 'East Outside';
L = legend('1920s','1930s','1940s','1950s','1960s','1970s','1980s','1990s','2000s','2010s','2020s',...
    '2030s','2040s','2050s','2060s','2070s'); L.Location = 'East Outside';
% for i = 1:84
% p = plot(1:7,-mmonthly_iceprod_20C(i,:)); p.Color = cmap1t(i,:);
% hold on;
% end
% for i = 1:73
% p = plot(1:7,-mmonthly_iceprod_RCP85(i,:)); p.Color = cmap1t(i+84,:);
% end
grid on;
title('Sea ice production through the winter season')
ax1.FontName = 'Helvetica'; ax1.FontSize = 13;
ax2 = axes('Position',[0.3375,0.18,0.3165,0.3*0.8]);
yyaxis left
p = plot(yrs,mpeak_iceprod_20C); p.LineWidth = 2; p.LineStyle = '-'; hold on;
p = plot(yrs85,mpeak_iceprod_RCP85); p.LineWidth = 2; p.LineStyle = '-';
p = plot([2006 2006],[2 5]); p.LineStyle = '-'; p.Color = [0.7 0.7 0.7];
t1 = text(1922,4.65,'Peak timing'); t1.Color = [0    0.4470    0.7410]; t1.FontName = 'Helvetica'; t1.FontSize = 12; t1.FontWeight = 'bold';

ax2.YTick = 2:5; ax2.YLim = [2 5]; ax2.XLim = [1920 2080];
ax2.YTickLabel = {'Nov','Dec','Jan','Feb'};
yyaxis right
p = plot(yrs,mpeak_iceprodval_20C/1e9); p.LineWidth = 2; p.LineStyle = '-';
p = plot(yrs85,mpeak_iceprodval_RCP85/1e9); p.LineWidth = 2; p.LineStyle = '-';
ax2.YLim = [300 550]; ax2.YLabel.String = 'km^3/month'
t2 = text(2020,520,'Peak value'); t2.Color = [0.8500    0.3250    0.0980]; t2.FontName = 'Helvetica'; t2.FontSize = 12; t2.FontWeight = 'bold';
ax2.FontName = 'Helvetica'; ax2.FontSize = 12;


% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/IceProduction_wintercycles.fig';
% saveas(gcf,fn);
% 
set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/IceProduction_wintercycles_40.pdf 

%%
figure;
ax1 = axes('Position',[0.125,0.1,0.8,0.8]);
for i = 1:86
p = plot(1:7,-mmonthly_I2O_20C(i,:)); p.Color = cmap1t(i,:);
hold on;
end
for i = 1:73
p = plot(1:7,-mmonthly_I2O_RCP85(i,:)); p.Color = cmap1t(i+84,:);
end
title('Ice mass balance')
ax1.FontName = 'Helvetica'; ax1.FontSize = 13; 

figure;
ax1 = axes('Position',[0.125,0.1,0.8,0.8]);
for i = 1:86
p = plot(1:7,-mmonthly_iceprod_20C(i,:)); p.Color = cmap1t(i,:);
hold on;
end
for i = 1:73
p = plot(1:7,-mmonthly_iceprod_RCP85(i,:)); p.Color = cmap1t(i+84,:);
end
title('Ice production')
ax1.FontName = 'Helvetica'; ax1.FontSize = 13;

figure;
ax1 = axes('Position',[0.125,0.1,0.8,0.8]);
for i = 1:86
p = plot(1:7,-mmonthly_melt_20C(i,:)); p.Color = cmap1t(i,:);
hold on;
end
for i = 1:73
p = plot(1:7,-mmonthly_melt_RCP85(i,:)); p.Color = cmap1t(i+84,:);
end
title('Melt')
ax1.FontName = 'Helvetica'; ax1.FontSize = 13;