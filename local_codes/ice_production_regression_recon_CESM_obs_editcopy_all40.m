% edit copy ice production: CESM and real world recon
% try the following:
% 2. subsample the ensemble to get estimate of error on the regressions
% coefficients
% 3. Test ensemble mean reconstruction with individual members
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


load CESM_T10m_20C_RCP85_all40.mat
load CESM_iceprod_vars_40.mat
% load CESM_Tair_20C_RCP85.mat SAT_20C SAT_RCP85



total_area = 1.2543e12;
SIAdeficit_effect_20C = total_area - SepSIA_20C;
SIAdeficit_effect_RCP85 = total_area - SepSIA_RCP85;

Sep_10mT_20C = squeeze(monthly_T10m_20C(:,1,:));
Sep_10mT_RCP85 = squeeze(monthly_T10m_RCP85(:,1,:));

compdiv_20C = pos_areadiv_20C - net_areadiv_20C;  % this is the area of divergence that is compensated by convergence and does not contribute to net area expansion
compdiv_RCP85 = pos_areadiv_RCP85 - net_areadiv_RCP85;



%% create timeseries

miceprod_20C = mean(iceprod_20C,2);
miceprod_RCP85 = mean(iceprod_RCP85,2);

mSAT_20C = mean(SAT_20C,2);
mSAT_RCP85 = mean(SAT_RCP85,2);

mSep_10mT_20C = mean(Sep_10mT_20C,2);
mSep_10mT_RCP85 = mean(Sep_10mT_RCP85,2);

deltaT_20C = -(SAT_20C+1.8);
deltaT_RCP85 = -(SAT_RCP85+1.8);

snow_effect_20C = deltaT_20C./hs_20C;   %deltaT/snowthickness
snow_effect_RCP85 = deltaT_RCP85./hs_RCP85;  
msnow_effect_20C = mean(snow_effect_20C,2);
msnow_effect_RCP85 = mean(snow_effect_RCP85,2);

SIAdeficit_effect_20C = deltaT_20C.*SIAdeficit_effect_20C;
SIAdeficit_effect_RCP85 = deltaT_RCP85.*SIAdeficit_effect_RCP85;
mSIAdeficit_effect_20C = mean(SIAdeficit_effect_20C,2);
mSIAdeficit_effect_RCP85 = mean(SIAdeficit_effect_RCP85,2);

netdiv_effect_20C = deltaT_20C.*net_areadiv_20C;
netdiv_effect_RCP85 = deltaT_RCP85.*net_areadiv_RCP85;
mnetdiv_effect_20C = mean(netdiv_effect_20C,2);
mnetdiv_effect_RCP85 = mean(netdiv_effect_RCP85,2);

totaldiv_effect_20C = deltaT_20C.*pos_areadiv_20C;    % using mean pos div here now rather than total div normalised by frz area days. May need to use mean pos div on freezing days however
totaldiv_effect_RCP85 = deltaT_RCP85.*pos_areadiv_RCP85;
mtotaldiv_effect_20C = mean(totaldiv_effect_20C,2);
mtotaldiv_effect_RCP85 = mean(totaldiv_effect_RCP85,2);


compdiv_effect_20C = deltaT_20C.*compdiv_20C;
compdiv_effect_RCP85 = deltaT_RCP85.*compdiv_RCP85;
mcompdiv_effect_20C = mean(compdiv_effect_20C,2);
mcompdiv_effect_RCP85 = mean(compdiv_effect_RCP85,2);

mFw_20C = mean(Fw_20C,2);
mFw_RCP85 = mean(Fw_RCP85,2);


% create full-length timeseries
snow_effect_full = cat(1,snow_effect_20C,snow_effect_RCP85);
msnow_effect_full = mean(snow_effect_full,2);

SIAdeficit_effect_full = cat(1,SIAdeficit_effect_20C,SIAdeficit_effect_RCP85);
mSIAdeficit_effect_full = mean(SIAdeficit_effect_full,2);

netdiv_effect_full = cat(1,netdiv_effect_20C,netdiv_effect_RCP85);
mnetdiv_effect_full = mean(netdiv_effect_full,2);

totaldiv_effect_full = cat(1,totaldiv_effect_20C,totaldiv_effect_RCP85);
mtotaldiv_effect_full = mean(totaldiv_effect_full,2);

compdiv_effect_full = cat(1,compdiv_effect_20C,compdiv_effect_RCP85);
mcompdiv_effect_full = mean(compdiv_effect_full,2);

Sep_10mT_full = cat(1,Sep_10mT_20C,Sep_10mT_RCP85);
mSep_10mT_full = mean(Sep_10mT_full,2);

SAT_full = cat(1,SAT_20C,SAT_RCP85);
mSAT_full = mean(SAT_full,2);

iceprod_full = cat(1,iceprod_20C,iceprod_RCP85);
miceprod_full = mean(iceprod_full,2);

Fw_full = cat(1,Fw_20C,Fw_RCP85);
mFw_full = mean(Fw_full,2);

%% detrended timeseries

for i = 1:40
diceprod_20C(:,i) = iceprod_20C(:,i) - miceprod_20C;
diceprod_RCP85(:,i) = iceprod_RCP85(:,i) - miceprod_RCP85;
diceprod_full(:,i) = iceprod_full(:,i) - miceprod_full;
dSAT_20C(:,i) = SAT_20C(:,i) - mSAT_20C;
dSAT_RCP85(:,i) = SAT_RCP85(:,i) - mSAT_RCP85;
dSAT_full(:,i) = SAT_full(:,i) - mSAT_full;
dSIAdeficit_effect_20C(:,i) = SIAdeficit_effect_20C(:,i) - mSIAdeficit_effect_20C;
dSIAdeficit_effect_RCP85(:,i) = SIAdeficit_effect_RCP85(:,i) - mSIAdeficit_effect_RCP85;
dSIAdeficit_effect_full(:,i) = SIAdeficit_effect_full(:,i) - mSIAdeficit_effect_full;
dnetdiv_effect_20C(:,i) = netdiv_effect_20C(:,i) - mnetdiv_effect_20C;
dnetdiv_effect_RCP85(:,i) = netdiv_effect_RCP85(:,i) - mnetdiv_effect_RCP85;
dnetdiv_effect_full(:,i) = netdiv_effect_full(:,i) - mnetdiv_effect_full;
dtotaldiv_effect_20C(:,i) = totaldiv_effect_20C(:,i) - mtotaldiv_effect_20C;
dtotaldiv_effect_RCP85(:,i) = totaldiv_effect_RCP85(:,i) - mtotaldiv_effect_RCP85;
dtotaldiv_effect_full(:,i) = totaldiv_effect_full(:,i) - mtotaldiv_effect_full;
dcompdiv_effect_20C(:,i) = compdiv_effect_20C(:,i) - mcompdiv_effect_20C;
dcompdiv_effect_RCP85(:,i) = compdiv_effect_RCP85(:,i) - mcompdiv_effect_RCP85;
dcompdiv_effect_full(:,i) = compdiv_effect_full(:,i) - mcompdiv_effect_full;
dsnow_effect_20C(:,i) = snow_effect_20C(:,i) - msnow_effect_20C;
dsnow_effect_RCP85(:,i) = snow_effect_RCP85(:,i) - msnow_effect_RCP85;
dsnow_effect_full(:,i) = snow_effect_full(:,i) - msnow_effect_full;
dSep_10mT_20C(:,i) = Sep_10mT_20C(:,i) - mSep_10mT_20C;
dSep_10mT_RCP85(:,i) = Sep_10mT_RCP85(:,i) - mSep_10mT_RCP85;
dSep_10mT_full(:,i) = Sep_10mT_full(:,i) - mSep_10mT_full;
dFw_20C(:,i) = Fw_20C(:,i) - mFw_20C;
dFw_RCP85(:,i) = Fw_RCP85(:,i) - mFw_RCP85;
dFw_full(:,i) = Fw_full(:,i) - mFw_full;
end

%% 1 sigma multi regression 20C...

% try normalising to get relative sensitivities
XSIAdeficit_effect = dSIAdeficit_effect_20C(:)/std(dSIAdeficit_effect_20C(:));
Xnetdiv_effect = dnetdiv_effect_20C(:)/std(dnetdiv_effect_20C(:));
Xcompdiv_effect = dcompdiv_effect_20C(:)/std(dcompdiv_effect_20C(:));
Xsnow_effect = dsnow_effect_20C(:)/std(dsnow_effect_20C(:));
XSep_10mT_effect = dSep_10mT_20C(:)/std(dSep_10mT_20C(:));
X_Fw_effect = dFw_20C(:)/std(dFw_20C(:));
% XFw = dFw_20C(:)/std(dFw_20C(:));
XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect,X_Fw_effect]; 
YIP = -diceprod_20C(:)/1e12;

lm_IP_20C = fitlm(XIP,YIP);

%% 1 sigma multi regression RCP85...

% try normalising to get relative sensitivities
XSIAdeficit_effect = dSIAdeficit_effect_RCP85(:)/std(dSIAdeficit_effect_20C(:));
Xnetdiv_effect = dnetdiv_effect_RCP85(:)/std(dnetdiv_effect_20C(:));
Xcompdiv_effect = dcompdiv_effect_RCP85(:)/std(dcompdiv_effect_20C(:));
Xsnow_effect = dsnow_effect_RCP85(:)/std(dsnow_effect_20C(:));
XSep_10mT_effect = dSep_10mT_RCP85(:)/std(dSep_10mT_20C(:));
X_Fw_effect = dFw_RCP85(:)/std(dFw_20C(:));

% XFw = dFw_RCP85(:)/std(dFw_RCP85(:));
XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect,X_Fw_effect]; 
YIP = -diceprod_RCP85(:)/1e12;

lm_IP_RCP85 = fitlm(XIP,YIP);

%% 1 sigma multi regression full...

% try normalising to get relative sensitivities
XSIAdeficit_effect = dSIAdeficit_effect_full(:)/std(dSIAdeficit_effect_20C(:));
Xnetdiv_effect = dnetdiv_effect_full(:)/std(dnetdiv_effect_20C(:));
Xcompdiv_effect = dcompdiv_effect_full(:)/std(dcompdiv_effect_20C(:));
Xsnow_effect = dsnow_effect_full(:)/std(dsnow_effect_20C(:));
XSep_10mT_effect = dSep_10mT_full(:)/std(dSep_10mT_20C(:));
X_Fw_effect = dFw_full(:)/std(dFw_20C(:)); 

% XFw = dFw_full(:)/std(dFw_full(:));
XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect,X_Fw_effect]; 
YIP = -diceprod_full(:)/1e12;

lm_IP_full = fitlm(XIP,YIP);

%% Figure: just the 1 sigma coefficients
c1 = [0 119 187]/256;
c2 = [51 187 238]/256;
c3 = [0 153 136]/256;
c4 = [238 119 51]/256;
c5 = [204 51 17]/256;
c6 = [238 51 119]/256;
c7 = [187 187 187]/256;
cIP = [c2;c3;c4;c7;c6];    % group the colours into arrays for each component
figure;
b_IP_20C = lm_IP_20C.Coefficients.Estimate(2:end);
b_IP_RCP85 = lm_IP_RCP85.Coefficients.Estimate(2:end);
b_IP_full = lm_IP_full.Coefficients.Estimate(2:end);
pos = [1,2,2.95;1,2,3.05;1,2,3;1,2,3;1,2,3];
ax1 = axes('Position',[0.125,0.1,0.775,0.8]);
name = {'(1/Snow depth) * \DeltaT','SIA deficit * \DeltaT','Net divergence * \DeltaT','Compensated divergence * \DeltaT','-Sep 10m T'};
for i = 1:length(b_IP_20C)
    p = plot(pos(i,:),[b_IP_20C(i),b_IP_RCP85(i),b_IP_full(i)]*1e3,'s','LineWidth',2,'MarkerSize',10); hold on;
    if i == 1
        txt = text(1.1,b_IP_20C(i)*1e3-0.01,char(name{i})); txt.Color = cIP(i,:); txt.FontName = 'Helvetica'; txt.FontSize = 12;
    elseif i==2
    txt = text(1.1,b_IP_20C(i)*1e3+0.01,char(name{i})); txt.Color = cIP(i,:); txt.FontName = 'Helvetica'; txt.FontSize = 12;
    else
    txt = text(1.1,b_IP_20C(i)*1e3,char(name{i})); txt.Color = cIP(i,:); txt.FontName = 'Helvetica'; txt.FontSize = 12;
    end
    p.Color = cIP(i,:); p.MarkerFaceColor = cIP(i,:);
end
xlim([0.5 3.5]); ax1.XTick = [1 2 3]; ax1.XTickLabel = {'20C','RCP85','Full'}; grid on; ylim([0 110]);
txt = text(0.82,5,strcat('R^2 = ',{' '},num2str(round(lm_IP_20C.Rsquared.Ordinary,2,'significant')))); txt.FontName = 'Helvetica'; txt.FontSize = 12; txt.Color = c5;
txt = text(1.82,5,strcat('R^2 = ',{' '},num2str(round(lm_IP_RCP85.Rsquared.Ordinary,2,'significant')))); txt.FontName = 'Helvetica'; txt.FontSize = 12; txt.Color = c5;
txt = text(2.82,5,strcat('R^2 = ',{' '},num2str(round(lm_IP_full.Rsquared.Ordinary,2,'significant')))); txt.FontName = 'Helvetica'; txt.FontSize = 12; txt.Color = c5;

t = title('Ice production coefficients, rescaled by 20C 1\sigma values'); ylabel('km^3'); 
ax1.FontName = 'Helvetica'; ax1.FontSize = 13;

 set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/IceProd_coefficients.pdf
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/IceProd_coefficients.png -m3

%% Plot in regression coefficients as bars
c_blue = [0.6 0.8 0.95];
c_orange = [0.95 0.8 0.6];
c_bar = [c_blue;c_orange;c1];
    cA = [0.6 0.6 1];    
name = {' (1/snow depth)*\DeltaT','Sep open water*\DeltaT','Net divergence*\DeltaT','Compensated\newline divergence*\DeltaT','Sep 10m T'};

figure;
ax = axes('Position',[0.125,0.125,0.8,0.8]);
% X = reshape(([0.9, 1, 1.1] + [0:length(b_IP_20C)-1])',1,length(b_IP_20C)*3); % creating an X axis for the bar chart
X = [0.9, 1, 1.1] + [0:length(b_IP_20C)]'; % creating an X axis for the bar chart
B = [b_IP_20C,b_IP_RCP85,b_IP_full]*1e3; % Putting the beta values into a matrix
E = [lm_IP_20C.Coefficients.SE(2:end),lm_IP_RCP85.Coefficients.SE(2:end),lm_IP_full.Coefficients.SE(2:end)]*1e3;
TE = [sum(lm_IP_20C.Coefficients.SE),sum(lm_IP_RCP85.Coefficients.SE),sum(lm_IP_full.Coefficients.SE)]*1e3;
tY = [-0.01,-0.01,-0.01,-0.01,0.01]*1e3;
for i = 1:length(b_IP_20C)
    for j = 1:3
p = bar(X(i,j),B(i,j),'ShowBaseline','off'); hold on;
er = errorbar(X(i,j),B(i,j),E(i,j),E(i,j)); er.Color = [0 0 0]; er.LineStyle = 'none';
p.FaceColor = c_bar(j,:); p.EdgeAlpha = 0;
p.FaceAlpha = cA(j); p.BarWidth = 0.08;
  
    end
    t = text(X(i,2),tY(i),name{i},'HorizontalAlignment','Center','Rotation',30,'FontSize',12);
end
L = legend('20C','RCP8.5','Full'); L.Location = 'NorthEast';
box off; grid on;
ax.XTick = [];
ax.FontName = 'Helvetica'; ax.FontSize = 13;
t = title('Ice production coefficients, rescaled by 20C 1 \sigma values'); ylabel('km^3'); 
% print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/IceProd_coefficients_bars.pdf','-dpdf','-bestfit')

%% detrended timeseries

for i = 1:40
diceprod_20C_mm(:,i) = iceprod_20C(:,i) - movmean(miceprod_20C,10);
diceprod_RCP85_mm(:,i) = iceprod_RCP85(:,i) - movmean(miceprod_RCP85,10);
diceprod_full_mm(:,i) = iceprod_full(:,i) - movmean(miceprod_full,10);
dSAT_20C_mm(:,i) = SAT_20C(:,i) - movmean(mSAT_20C,10);
dSAT_RCP85_mm(:,i) = SAT_RCP85(:,i) - movmean(mSAT_RCP85,10);
dSAT_full_mm(:,i) = SAT_full(:,i) - movmean(mSAT_full,10);
dSIAdeficit_effect_20C_mm(:,i) = SIAdeficit_effect_20C(:,i) - movmean(mSIAdeficit_effect_20C,10);
dSIAdeficit_effect_RCP85_mm(:,i) = SIAdeficit_effect_RCP85(:,i) - movmean(mSIAdeficit_effect_RCP85,10);
dSIAdeficit_effect_full_mm(:,i) = SIAdeficit_effect_full(:,i) - movmean(mSIAdeficit_effect_full,10);
dnetdiv_effect_20C_mm(:,i) = netdiv_effect_20C(:,i) - movmean(mnetdiv_effect_20C,10);
dnetdiv_effect_RCP85_mm(:,i) = netdiv_effect_RCP85(:,i) - movmean(mnetdiv_effect_RCP85,10);
dnetdiv_effect_full_mm(:,i) = netdiv_effect_full(:,i) - movmean(mnetdiv_effect_full,10);
dtotaldiv_effect_20C_mm(:,i) = totaldiv_effect_20C(:,i) - movmean(mtotaldiv_effect_20C,10);
dtotaldiv_effect_RCP85_mm(:,i) = totaldiv_effect_RCP85(:,i) - movmean(mtotaldiv_effect_RCP85,10);
dtotaldiv_effect_full_mm(:,i) = totaldiv_effect_full(:,i) - movmean(mtotaldiv_effect_full,10);
dcompdiv_effect_20C_mm(:,i) = compdiv_effect_20C(:,i) - movmean(mcompdiv_effect_20C,10);
dcompdiv_effect_RCP85_mm(:,i) = compdiv_effect_RCP85(:,i) - movmean(mcompdiv_effect_RCP85,10);
dcompdiv_effect_full_mm(:,i) = compdiv_effect_full(:,i) - movmean(mcompdiv_effect_full,10);
dsnow_effect_20C_mm(:,i) = snow_effect_20C(:,i) - movmean(msnow_effect_20C,10);
dsnow_effect_RCP85_mm(:,i) = snow_effect_RCP85(:,i) - movmean(msnow_effect_RCP85,10);
dsnow_effect_full_mm(:,i) = snow_effect_full(:,i) - movmean(msnow_effect_full,10);
dSep_10mT_20C_mm(:,i) = Sep_10mT_20C(:,i) - movmean(mSep_10mT_20C,10);
dSep_10mT_RCP85_mm(:,i) = Sep_10mT_RCP85(:,i) - movmean(mSep_10mT_RCP85,10);
dSep_10mT_full_mm(:,i) = Sep_10mT_full(:,i) - movmean(mSep_10mT_full,10);
dFw_20C_mm(:,i) = Fw_20C(:,i) - movmean(mFw_20C,10);
dFw_RCP85_mm(:,i) = Fw_RCP85(:,i) - movmean(mFw_RCP85,10);
dFw_full_mm(:,i) = Fw_full(:,i) - movmean(mFw_full,10);
end
%% 1 sigma multi regression 20C... FOR EACH ENS MEMBER

% NOTE: reintroduce _mm to use the deviations from the moving ensemble mean
% rather than just deviations from the ensemble mean
for i = 1:40
% try normalising to get relative sensitivities
XSIAdeficit_effect = dSIAdeficit_effect_20C(:,i)/std(dSIAdeficit_effect_20C(:,i));
Xnetdiv_effect = dnetdiv_effect_20C(:,i)/std(dnetdiv_effect_20C(:,i));
Xcompdiv_effect = dcompdiv_effect_20C(:,i)/std(dcompdiv_effect_20C(:,i));
Xsnow_effect = dsnow_effect_20C(:,i)/std(dsnow_effect_20C(:,i));
XSep_10mT_effect = dSep_10mT_20C(:,i)/std(dSep_10mT_20C(:,i));
X_Fw_effect = dFw_20C(:,i)/std(dFw_20C(:,i));
% XFw = dFw_20C(:)/std(dFw_20C(:));
XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect]; 
YIP = -diceprod_20C(:,i)/1e12;

lm_IP_20C = fitlm(XIP,YIP);
b_ens_IP_20C(:,i) = lm_IP_20C.Coefficients.Estimate(2:end);
err_ens_IP_20C(:,i) = lm_IP_20C.Coefficients.SE(2:end);
end

%% 1 sigma multi regression RCP85...
for i = 1:40
% try normalising to get relative sensitivities
XSIAdeficit_effect = dSIAdeficit_effect_RCP85(:,i)/std(dSIAdeficit_effect_20C(:,i));
Xnetdiv_effect = dnetdiv_effect_RCP85(:,i)/std(dnetdiv_effect_20C(:,i));
Xcompdiv_effect = dcompdiv_effect_RCP85(:,i)/std(dcompdiv_effect_20C(:,i));
Xsnow_effect = dsnow_effect_RCP85(:,i)/std(dsnow_effect_20C(:,i));
XSep_10mT_effect = dSep_10mT_RCP85(:,i)/std(dSep_10mT_20C(:,i));
X_Fw_effect = dFw_RCP85(:,i)/std(dFw_20C(:,i));

% XFw = dFw_RCP85(:)/std(dFw_RCP85(:));
XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect]; 
YIP = -diceprod_RCP85(:,i)/1e12;

lm_IP_RCP85 = fitlm(XIP,YIP);
b_ens_IP_RCP85(:,i) = lm_IP_RCP85.Coefficients.Estimate(2:end);
err_ens_IP_RCP85(:,i) = lm_IP_RCP85.Coefficients.SE(2:end);
end

%% 1 sigma multi regression full...
for i = 1:40
% try normalising to get relative sensitivities
XSIAdeficit_effect = dSIAdeficit_effect_full(:,i)/std(dSIAdeficit_effect_20C(:,i));
Xnetdiv_effect = dnetdiv_effect_full(:,i)/std(dnetdiv_effect_20C(:,i));
Xcompdiv_effect = dcompdiv_effect_full(:,i)/std(dcompdiv_effect_20C(:,i));
Xsnow_effect = dsnow_effect_full(:,i)/std(dsnow_effect_20C(:,i));
XSep_10mT_effect = dSep_10mT_full(:,i)/std(dSep_10mT_20C(:,i));
X_Fw_effect = dFw_full(:,i)/std(dFw_20C(:,i)); 

% XFw = dFw_full(:)/std(dFw_full(:));
XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect]; 
YIP = -diceprod_full(:,i)/1e12;

lm_IP_full = fitlm(XIP,YIP);
b_ens_IP_full(:,i) = lm_IP_full.Coefficients.Estimate(2:end);
err_ens_IP_full(:,i) = lm_IP_full.Coefficients.SE(2:end);
end

%% Plot in regression coefficients as bars with uncertainty estimate
c_blue = [0.6 0.8 0.95];
c_orange = [0.95 0.8 0.6];
c_bar = [c_blue;c_orange;c1];
    cA = [0.6 0.6 1];    
name = {' (1/snow depth)*\DeltaT','Sep open water*\DeltaT','Net divergence*\DeltaT','Compensated\newline divergence*\DeltaT','Sep 10m T'};

figure;
ax = axes('Position',[0.125,0.15,0.8,0.775]);
% X = reshape(([0.9, 1, 1.1] + [0:length(b_IP_20C)-1])',1,length(b_IP_20C)*3); % creating an X axis for the bar chart
X = [0.9, 1, 1.1] + [0:length(b_IP_20C)]'; % creating an X axis for the bar chart
B = [b_IP_20C,b_IP_RCP85,b_IP_full]*1e3; % Putting the beta values into a matrix
tY = [-0.0125,-0.0125,-0.0125,-0.0125,0.01]*1e3;
b_ens = cat(3,b_ens_IP_20C,b_ens_IP_RCP85,b_ens_IP_full)*1e3;
for i = 1:length(b_IP_20C)
    for j = 1:3
p = bar(X(i,j),B(i,j),'ShowBaseline','off');
p.FaceColor = c_bar(j,:); p.EdgeAlpha = 0;
p.FaceAlpha = cA(j); p.BarWidth = 0.08;
    hold on;
    end
%     t = text(X(i,2),tY(i),name{i},'HorizontalAlignment','Center','Rotation',30,'FontSize',12)
end
L = legend('20C','RCP8.5','Full'); L.Location = 'NorthEast';
for i = 1:length(b_IP_20C); for j = 1:3
pd = plot(X(i,j)*ones(1,40),b_ens(i,:,j),'Marker','.','LineStyle','none','Color',c7);
end; end
% find the mean and std of the individual ensemble estimates
b_ens_mean = squeeze(mean(b_ens,2));
b_ens_std = squeeze(std(b_ens,[],2));
for i = 1:length(b_IP_20C); for j = 1:3
errpl = plot(X(i,j),b_ens_mean(i,j),'Marker','.','LineStyle','none','Color','k');
er = errorbar(X(i,j),b_ens_mean(i,j),b_ens_std(i,j),b_ens_std(i,j)); er.Color = [0 0 0]; er.LineStyle = 'none';
end; end

p = plot([0.5,length(b_IP_20C)+0.5],[0, 0],'color',[0.5 0.5 0.5]);
box off; grid on;

% ax.XTick = [];
ax.XTick = [1:length(b_IP_20C)];
ax.XTickLabel = {'(1/h_s)*\DeltaT','Sep open\newline water*\DeltaT','   Net\newline div.*\DeltaT',' Comp.\newline div.*\DeltaT','Sep 10m T'};
% ax.XTickLabelRotation = 30;
ax.FontName = 'Helvetica'; ax.FontSize = 13;
t = title('Ice production coefficients, rescaled by 20C 1 \sigma values'); ylabel('km^3'); 
% print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/IceProd_coefficients_bars_uncertainty_all40.pdf','-dpdf','-bestfit')

%%  plot reconstruction with separate window for individual regressors

XSIAdeficit_effect = dSIAdeficit_effect_full(:)/1e12;
Xnetdiv_effect = dnetdiv_effect_full(:)/1e12;
Xcompdiv_effect = dcompdiv_effect_full(:)/1e12;
Xsnow_effect = dsnow_effect_full(:);
XSep_10mT = dSep_10mT_full(:);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT]; 
YIP = -diceprod_full(:)/1e12;

lm_mrIP_full = fitlm(XIP,YIP); 
beta = lm_mrIP_full.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

iceprod_est_full_byfull = -mean(iceprod_full(:)) + beta(1)...
    + beta(2).*(msnow_effect_full-mean(msnow_effect_full))...
    + beta(3).*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))...
    + beta(4)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))...
    + beta(5).*(mcompdiv_effect_full-mean(mcompdiv_effect_full))...
    + beta(6)*(mSep_10mT_full-mean(mSep_10mT_full));
clearvars beta
XSIAdeficit_effect = dSIAdeficit_effect_20C(:)/1e12;
Xnetdiv_effect = dnetdiv_effect_20C(:)/1e12;
Xcompdiv_effect = dcompdiv_effect_20C(:)/1e12;
Xsnow_effect = dsnow_effect_20C(:);
XSep_10mT = dSep_10mT_20C(:);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT]; 
YIP = -diceprod_20C(:)/1e12;

lm_mrIP_20C = fitlm(XIP,YIP); 
beta = lm_mrIP_20C.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

iceprod_est_full_by20C = -mean(iceprod_full(:)) + beta(1)...
    + beta(2).*(msnow_effect_full-mean(msnow_effect_full))...
    + beta(3).*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))...
    + beta(4)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))...
    + beta(5).*(mcompdiv_effect_full-mean(mcompdiv_effect_full))...
    + beta(6)*(mSep_10mT_full-mean(mSep_10mT_full));
clearvars beta
XSIAdeficit_effect = dSIAdeficit_effect_RCP85(:)/1e12;
Xnetdiv_effect = dnetdiv_effect_RCP85(:)/1e12;
Xcompdiv_effect = dcompdiv_effect_RCP85(:)/1e12;
Xsnow_effect = dsnow_effect_RCP85(:);
XSep_10mT = dSep_10mT_RCP85(:);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT]; 
YIP = -diceprod_RCP85(:)/1e12;

lm_mrIP_RCP85 = fitlm(XIP,YIP); 
beta = lm_mrIP_RCP85.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

iceprod_est_full_byRCP85 = -mean(iceprod_full(:)) + beta(1)...
    + beta(2).*(msnow_effect_full-mean(msnow_effect_full))...
    + beta(3).*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))...
    + beta(4)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))...
    + beta(5).*(mcompdiv_effect_full-mean(mcompdiv_effect_full))...
    + beta(6)*(mSep_10mT_full-mean(mSep_10mT_full));
clearvars beta
c_bluet = [0.6 0.8 0.95 0.6];
c_oranget = [0.95 0.8 0.6 0.6];
% plot
figure;
ax1 = axes('Position',[0.125,0.35,0.75,0.6]);
p = plot(yrs_full, 1e-12*((-1*miceprod_full)),'color','k'); p.LineWidth = 2;  hold on;
p = plot(yrs_full, 1e-12*(iceprod_est_full_byfull-1*mean(miceprod_full(:))-mean(iceprod_est_full_byfull(:))),'color',c1); p.LineWidth = 2; 
p = plot(yrs_full, 1e-12*(iceprod_est_full_by20C-1*mean(miceprod_full(:))-mean(iceprod_est_full_by20C(:))),'color',c_bluet); p.LineWidth = 2; 
p = plot(yrs_full, 1e-12*(iceprod_est_full_byRCP85-1*mean(miceprod_full(:))-mean(iceprod_est_full_byRCP85(:))),'color',c_oranget); p.LineWidth = 2; 
L = legend('ensemble mean','est: full record','est: 20C data','est: RCP8.5 data'); L.Location = 'SouthWest';
p = plot([2006.5 2006.5],[1.4 2.4]); p.LineStyle = '-'; p.Color = [0.3 0.3 0.3];
tx1 = text(1990,1.95,'20C'); tx2 = text(2011,1.95,'RCP8.5'); tx1.FontName = 'Helvetica'; tx2.FontName = 'Helvetica'; tx1.FontSize = 13; tx2.FontSize = 13; tx1.Color = [0.3 0.3 0.3]; tx2.Color = [0.3 0.3 0.3];
grid on;
t = title('Linear model for winter ice production'); t.FontName = 'Helvetica'; t.FontWeight = 'bold';
 
ax1.FontSize = 12; ax1.FontName = 'Helvetica'; ax1.XTickLabel = [];
ax1.YLim = [1.4 2.4];


ax3 = axes('Position',[0.125,0.075,0.75,0.24]);
beta = lm_mrIP_full.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

p = plot(yrs_full, 1e-12*(beta(2)*(msnow_effect_full-mean(msnow_effect_full))-beta(2)*(msnow_effect_full(1)-mean(msnow_effect_full)))); p.Color = c2; hold on; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(3)*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))-beta(3)*(mSIAdeficit_effect_full(1)-mean(mSIAdeficit_effect_full)))); p.Color = c3; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(4)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))-beta(4)*(mnetdiv_effect_full(1)-mean(mnetdiv_effect_full)))); p.Color = c4; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(5)*(mcompdiv_effect_full-mean(mcompdiv_effect_full))-beta(5)*(mcompdiv_effect_full(1)-mean(mcompdiv_effect_full)))); p.Color = c7; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(6)*(mSep_10mT_full-mean(mSep_10mT_full))-beta(6)*(mSep_10mT_full(1)-mean(mSep_10mT_full)))); p.Color = c6; p.LineWidth = 1.5;
L1 = legend('(1/Snow depth) *\DeltaT','Sep open water *\DeltaT','Net divergence *\DeltaT','Comp. divergence *\DeltaT','Sep 10 m T'); L1.Position = L.Position + [0.35,0.01,L1.Position(3)-L.Position(3),L1.Position(4)-L.Position(4)];
p = plot([2006.5 2006.5],[-0.2 0.2]); p.LineStyle = '-'; p.Color = [0.3 0.3 0.3];
ax3.YLim = [-0.2 0.2];
grid on;
% t = title('Winter ice production linear model trained on full record'); t.FontName = 'Futura'; t.FontWeight = 'bold';
% ylabel('1000 km^3');
ax3.FontSize = 12; ax3.FontName = 'Helvetica';

sy = suplabel('Ice production, 1000 km^3','y'); sy.FontSize = 12; sy.FontName = 'Helvetica';

t1 = text(ax1,1925,2.325,'Full reconstruction'); t1.FontName = 'Helvetica'; t1.FontSize = 13; t1.FontWeight = 'bold';
t2 = text(ax3,1925,0.15,'Constituent parts'); t2.FontName = 'Helvetica'; t2.FontSize = 13; t2.FontWeight = 'bold';
t3 = text(ax1,1910,2.3,'a)'); t3.FontName = 'Helvetica'; t3.FontSize = 13; t3.FontWeight = 'bold';
t4 = text(ax3,1910,0.125,'b)'); t4.FontName = 'Helvetica'; t4.FontSize = 13; t4.FontWeight = 'bold';

set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/Ice_production_recon_2_panelsv2.pdf

% print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/Ice_production_recon_2_panelsv2_edit_40.pdf','-dpdf','-bestfit');

%% Reconstruction with individual ensemble member estimates as uncertainty window
clearvars beta_ens_IP_20C beta_ens_IP_RCP85 beta_ens_IP_full
for i = 1:40
% 20C
XSIAdeficit_effect = dSIAdeficit_effect_20C(:,i)/1e12;
Xnetdiv_effect = dnetdiv_effect_20C(:,i)/1e12;
Xcompdiv_effect = dcompdiv_effect_20C(:,i)/1e12;
Xsnow_effect = dsnow_effect_20C(:,i);
XSep_10mT_effect = dSep_10mT_20C(:,i);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect]; 
YIP = -diceprod_20C(:,i)/1e12;
% -mean(iceprod_full(:))
lm_IP_20C = fitlm(XIP,YIP);
beta_ens_IP_20C(:,i) = lm_IP_20C.Coefficients.Estimate; beta_ens_IP_20C(3:5,i) = beta_ens_IP_20C(3:5,i)/1e12;  beta_ens_IP_20C(:,i) = beta_ens_IP_20C(:,i)*1e12;
iceprod_est_full_by20C_allensmem(:,i) = beta_ens_IP_20C(1,i)...
    + beta_ens_IP_20C(2,i).*(msnow_effect_full-mean(msnow_effect_full))...
    + beta_ens_IP_20C(3,i).*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))...
    + beta_ens_IP_20C(4,i)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))...
    + beta_ens_IP_20C(5,i).*(mcompdiv_effect_full-mean(mcompdiv_effect_full))...
    + beta_ens_IP_20C(6,i)*(mSep_10mT_full-mean(mSep_10mT_full));

%RCP85
XSIAdeficit_effect = dSIAdeficit_effect_RCP85(:,i)/1e12;
Xnetdiv_effect = dnetdiv_effect_RCP85(:,i)/1e12;
Xcompdiv_effect = dcompdiv_effect_RCP85(:,i)/1e12;
Xsnow_effect = dsnow_effect_RCP85(:,i);
XSep_10mT_effect = dSep_10mT_RCP85(:,i);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect]; 
YIP = -diceprod_RCP85(:,i)/1e12;

lm_IP_RCP85 = fitlm(XIP,YIP);
beta_ens_IP_RCP85(:,i) = lm_IP_RCP85.Coefficients.Estimate; beta_ens_IP_RCP85(3:5,i) = beta_ens_IP_RCP85(3:5,i)/1e12;  beta_ens_IP_RCP85(:,i) = beta_ens_IP_RCP85(:,i)*1e12;
iceprod_est_full_byRCP85_allensmem(:,i) = beta_ens_IP_RCP85(1,i)...
    + beta_ens_IP_RCP85(2,i).*(msnow_effect_full-mean(msnow_effect_full))...
    + beta_ens_IP_RCP85(3,i).*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))...
    + beta_ens_IP_RCP85(4,i)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))...
    + beta_ens_IP_RCP85(5,i).*(mcompdiv_effect_full-mean(mcompdiv_effect_full))...
    + beta_ens_IP_RCP85(6,i)*(mSep_10mT_full-mean(mSep_10mT_full));

%full
XSIAdeficit_effect = dSIAdeficit_effect_full(:,i)/1e12;
Xnetdiv_effect = dnetdiv_effect_full(:,i)/1e12;
Xcompdiv_effect = dcompdiv_effect_full(:,i)/1e12;
Xsnow_effect = dsnow_effect_full(:,i);
XSep_10mT_effect = dSep_10mT_full(:,i);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT_effect]; 
YIP = -diceprod_full(:,i)/1e12;

lm_IP_full = fitlm(XIP,YIP);
beta_ens_IP_full(:,i) = lm_IP_full.Coefficients.Estimate; beta_ens_IP_full(3:5,i) = beta_ens_IP_full(3:5,i)/1e12;  beta_ens_IP_full(:,i) = beta_ens_IP_full(:,i)*1e12;
iceprod_est_full_byfull_allensmem(:,i) = beta_ens_IP_full(1,i)...
    + beta_ens_IP_full(2,i).*(msnow_effect_full-mean(msnow_effect_full))...
    + beta_ens_IP_full(3,i).*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))...
    + beta_ens_IP_full(4,i)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))...
    + beta_ens_IP_full(5,i).*(mcompdiv_effect_full-mean(mcompdiv_effect_full))...
    + beta_ens_IP_full(6,i)*(mSep_10mT_full-mean(mSep_10mT_full));
end

%  find mean and std through years
    em_iceprod_est_full_by20C = mean(iceprod_est_full_by20C_allensmem,2);
    em_iceprod_est_full_byRCP85 = mean(iceprod_est_full_byRCP85_allensmem,2);
    em_iceprod_est_full_byfull = mean(iceprod_est_full_byfull_allensmem,2);
    
    esd_iceprod_est_full_by20C = std(iceprod_est_full_by20C_allensmem,[],2);
    esd_iceprod_est_full_byRCP85 = std(iceprod_est_full_byRCP85_allensmem,[],2);
    esd_iceprod_est_full_byfull = std(iceprod_est_full_byfull_allensmem,[],2);

figure;
ax1 = axes('Position',[0.125,0.35,0.75,0.6]);
p = plot(yrs_full, 1e-12*((-1*miceprod_full)),'color','k'); p.LineWidth = 2;  hold on;
[h1, hp1] = boundedline(yrs_full,1e-12*(em_iceprod_est_full_byfull-1*mean(miceprod_full(:))-mean(em_iceprod_est_full_byfull(:))),1e-12*esd_iceprod_est_full_byfull); h1.Color = c1; hp1.FaceColor = c1; hp1.FaceAlpha = 0.5; 
[h2, hp2] = boundedline(yrs_full,1e-12*(em_iceprod_est_full_by20C-1*mean(miceprod_full(:))-mean(em_iceprod_est_full_by20C(:))),1e-12*esd_iceprod_est_full_by20C); h2.Color = c_blue; hp2.FaceColor = c_blue; hp2.FaceAlpha = 0.5; 
[h3, hp3] = boundedline(yrs_full,1e-12*(em_iceprod_est_full_byRCP85-1*mean(miceprod_full(:))-mean(em_iceprod_est_full_byRCP85(:))),1e-12*esd_iceprod_est_full_byRCP85); h3.Color = c_orange; hp3.FaceColor = c_orange; hp3.FaceAlpha = 0.5;

figure;
ax1 = axes('Position',[0.125,0.35,0.75,0.6]);

p = plot(yrs_full, 1e-12*((-1*miceprod_full)),'color','k'); p.LineWidth = 2;  hold on;
p = plot(yrs_full, 1e-12*(iceprod_est_full_byfull-1*mean(miceprod_full(:))-mean(iceprod_est_full_byfull(:))),'color',c1); p.LineWidth = 2; 
p = plot(yrs_full, 1e-12*(iceprod_est_full_by20C-1*mean(miceprod_full(:))-mean(iceprod_est_full_by20C(:))),'color',c_bluet); p.LineWidth = 2; 
p = plot(yrs_full, 1e-12*(iceprod_est_full_byRCP85-1*mean(miceprod_full(:))-mean(iceprod_est_full_byRCP85(:))),'color',c_oranget); p.LineWidth = 2; 
L = legend('ensemble mean','est: full record','est: 20C data','est: RCP8.5 data'); L.Location = 'SouthWest';
L.AutoUpdate = 'off';
[h1, hp1] = boundedline(yrs_full,1e-12*(em_iceprod_est_full_byfull-1*mean(miceprod_full(:))-mean(em_iceprod_est_full_byfull(:))),1e-12*esd_iceprod_est_full_byfull); h1.Color = c1; hp1.FaceColor = c1; hp1.FaceAlpha = 0.5; h1.LineStyle = 'none';
[h2, hp2] = boundedline(yrs_full,1e-12*(em_iceprod_est_full_by20C-1*mean(miceprod_full(:))-mean(em_iceprod_est_full_by20C(:))),1e-12*esd_iceprod_est_full_by20C); h2.Color = c_blue; hp2.FaceColor = c_blue; hp2.FaceAlpha = 0.5; h2.LineStyle = 'none';
[h3, hp3] = boundedline(yrs_full,1e-12*(em_iceprod_est_full_byRCP85-1*mean(miceprod_full(:))-mean(em_iceprod_est_full_byRCP85(:))),1e-12*esd_iceprod_est_full_byRCP85); h3.Color = c_orange; hp3.FaceColor = c_orange; hp3.FaceAlpha = 0.5; h3.LineStyle = 'none';
p = plot(yrs_full, 1e-12*((-1*miceprod_full)),'color','k'); p.LineWidth = 2;  hold on;
p = plot(yrs_full, 1e-12*(iceprod_est_full_byfull-1*mean(miceprod_full(:))-mean(iceprod_est_full_byfull(:))),'color',c1); p.LineWidth = 2; 
p = plot(yrs_full, 1e-12*(iceprod_est_full_by20C-1*mean(miceprod_full(:))-mean(iceprod_est_full_by20C(:))),'color',c_bluet); p.LineWidth = 2; 
p = plot(yrs_full, 1e-12*(iceprod_est_full_byRCP85-1*mean(miceprod_full(:))-mean(iceprod_est_full_byRCP85(:))),'color',c_oranget); p.LineWidth = 2;

p = plot([2006.5 2006.5],[1.4 2.4]); p.LineStyle = '-'; p.Color = [0.3 0.3 0.3];

tx1 = text(1990,1.95,'20C'); tx2 = text(2011,1.95,'RCP8.5'); tx1.FontName = 'Helvetica'; tx2.FontName = 'Helvetica'; tx1.FontSize = 13; tx2.FontSize = 13; tx1.Color = [0.3 0.3 0.3]; tx2.Color = [0.3 0.3 0.3];
grid on;
t = title('Application of linear model to forced changes in ice production'); t.FontName = 'Helvetica'; t.FontWeight = 'bold';
 
ax1.FontSize = 12; ax1.FontName = 'Helvetica'; ax1.XTickLabel = [];
ax1.YLim = [1.4 2.4];


ax3 = axes('Position',[0.125,0.075,0.75,0.24]);
beta = lm_mrIP_full.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

p = plot(yrs_full, 1e-12*(beta(2)*(msnow_effect_full-mean(msnow_effect_full))-beta(2)*(msnow_effect_full(1)-mean(msnow_effect_full)))); p.Color = c2; hold on; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(3)*(mSIAdeficit_effect_full-mean(mSIAdeficit_effect_full))-beta(3)*(mSIAdeficit_effect_full(1)-mean(mSIAdeficit_effect_full)))); p.Color = c3; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(4)*(mnetdiv_effect_full-mean(mnetdiv_effect_full))-beta(4)*(mnetdiv_effect_full(1)-mean(mnetdiv_effect_full)))); p.Color = c4; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(5)*(mcompdiv_effect_full-mean(mcompdiv_effect_full))-beta(5)*(mcompdiv_effect_full(1)-mean(mcompdiv_effect_full)))); p.Color = c7; p.LineWidth = 1.5;
p = plot(yrs_full, 1e-12*(beta(6)*(mSep_10mT_full-mean(mSep_10mT_full))-beta(6)*(mSep_10mT_full(1)-mean(mSep_10mT_full)))); p.Color = c6; p.LineWidth = 1.5;
L1 = legend('\beta_1(1/Snow depth)\DeltaT','\beta_2Sep open water\DeltaT','\beta_3Net divergence\DeltaT','\beta_4Comp. divergence\DeltaT','\beta_5Sep SST'); L1.Position = L.Position + [0.345,-0.07,L1.Position(3)-L.Position(3),L1.Position(4)-L.Position(4)];
p = plot([2006.5 2006.5],[-0.2 0.2]); p.LineStyle = '-'; p.Color = [0.3 0.3 0.3];
ax3.YLim = [-0.2 0.2];
grid on;
% t = title('Winter ice production linear model trained on full record'); t.FontName = 'Futura'; t.FontWeight = 'bold';
% ylabel('1000 km^3');
ax3.FontSize = 12; ax3.FontName = 'Helvetica';

sy = suplabel('Winter ice production, 1000 km^3','y'); sy.FontSize = 12; sy.FontName = 'Helvetica';

t1 = text(ax1,1925,2.325,'Full reconstruction'); t1.FontName = 'Helvetica'; t1.FontSize = 13; t1.FontWeight = 'bold';
t2 = text(ax3,1925,0.15,'Constituent parts'); t2.FontName = 'Helvetica'; t2.FontSize = 13; t2.FontWeight = 'bold';
t3 = text(ax1,1910,2.3,'a)'); t3.FontName = 'Helvetica'; t3.FontSize = 13; t3.FontWeight = 'bold';
t4 = text(ax3,1910,0.125,'b)'); t4.FontName = 'Helvetica'; t4.FontSize = 13; t4.FontWeight = 'bold';

set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/Ice_production_recon_2_panelsv3.pdf

% print(gcf,'~/Documents/MATLAB/project/data/CESM/figures/Ice_production_recon_2_panelsv4_40.pdf','-dpdf','-bestfit');

%% Plot scatterplots alongside regression coefficients plot

% first: regression coefficients
figure;
ax = axes('Position',[0.1,0.15,0.55,0.775]);
% X = reshape(([0.9, 1, 1.1] + [0:length(b_IP_20C)-1])',1,length(b_IP_20C)*3); % creating an X axis for the bar chart
X = [0.9, 1, 1.1] + [0:length(b_IP_20C)]'; % creating an X axis for the bar chart

for i = 1:length(b_IP_20C)
    for j = 1:3
p = bar(X(i,j),B(i,j),'ShowBaseline','off');
p.FaceColor = c_bar(j,:); p.EdgeAlpha = 0;
p.FaceAlpha = cA(j); p.BarWidth = 0.08;
    hold on;
    end
%     t = text(X(i,2),tY(i),name{i},'HorizontalAlignment','Center','Rotation',30,'FontSize',12)
end
L = legend('20C','RCP8.5','Full'); L.Location = 'NorthEast';
for i = 1:length(b_IP_20C); for j = 1:3
pd = plot(X(i,j)*ones(1,40),b_ens(i,:,j),'Marker','.','LineStyle','none','Color',c7);
end; end
% mean and std of the individual ensemble estimates
for i = 1:length(b_IP_20C); for j = 1:3
errpl = plot(X(i,j),b_ens_mean(i,j),'Marker','.','LineStyle','none','Color','k');
er = errorbar(X(i,j),b_ens_mean(i,j),b_ens_std(i,j),b_ens_std(i,j)); er.Color = [0 0 0]; er.LineStyle = 'none';
end; end
p = plot([0.5,length(b_IP_20C)+0.5],[0, 0],'color',[0.5 0.5 0.5]);
box on; grid on;

% ax.XTick = [];
ax.XTick = [1:length(b_IP_20C)];
ax.XTickLabel = {'\beta_1\newline (1/h_s)*\DeltaT','\beta_2\newline Sep open\newline water*\DeltaT','\beta_3\newline   Net\newline div.*\DeltaT','\beta_4\newline Comp.\newline div.*\DeltaT','\beta_5\newline Sep 10m T'};
ax.XTickLabel = {'(1/h_s)*\DeltaT','Sep open\newline water*\DeltaT','   Net\newline div.*\DeltaT',' Comp.\newline div.*\DeltaT','Sep 10m T'};
ax.XTickLabel = {'\beta_1','\beta_2','\beta_3','\beta_4','\beta_5'};
% ax.XTickLabelRotation = 30;
ax.FontName = 'Helvetica'; ax.FontSize = 13;
t = title('Ice production coefficients, rescaled by 20C 1 \sigma values'); ylabel('Winter ice production, km^3'); t.Position = t.Position + [0 5 0];
ax.XLim = [0.5 5.5];


% 20C scatter plot
ax3 = axes('Position',[0.75,0.7,0.2,0.225]);

YIP = -diceprod_20C(:)/1e12;
beta = lm_mrIP_20C.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

iceprod_est_20C_alldata = beta(1)...
    + beta(2)*dsnow_effect_20C(:)...
    + beta(3)*dSIAdeficit_effect_20C(:)...
    + beta(4)*dnetdiv_effect_20C(:)...
    + beta(5)*dcompdiv_effect_20C(:)...
    + beta(6)*dSep_10mT_20C(:);

s1 = scatter(iceprod_est_20C_alldata/1e9,YIP*1e3); s1.MarkerEdgeColor = c_blue; s1.Marker = '.'; s1.MarkerEdgeAlpha = 0.4; s1.SizeData = 30;
hold on; grid on; xlim([-500 500]); ylim([-500 500]);
x = [-400:400];
lobf = polyfit(iceprod_est_20C_alldata/1e9,YIP*1e3,1);
p = plot(x,lobf(1)*x+lobf(2),'r');
t_lmf = title('Linear model fits'); t_lmf.Position = t_lmf.Position + [0 20 0];
ax3.XTickLabel = [];


% RCP85 scatter plot
ax2 = axes('Position',[0.75,0.425,0.2,0.225]);

YIP = -diceprod_RCP85(:)/1e12;
beta = lm_mrIP_RCP85.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

iceprod_est_RCP85_alldata = beta(1)...
    + beta(2)*dsnow_effect_RCP85(:)...
    + beta(3)*dSIAdeficit_effect_RCP85(:)...
    + beta(4)*dnetdiv_effect_RCP85(:)...
    + beta(5)*dcompdiv_effect_RCP85(:)...
    + beta(6)*dSep_10mT_RCP85(:);

s1 = scatter(iceprod_est_RCP85_alldata/1e9,YIP*1e3); s1.MarkerEdgeColor = c_orange; s1.Marker = '.'; s1.MarkerEdgeAlpha = 0.4; s1.SizeData = 30;
hold on; grid on; xlim([-500 500]); ylim([-500 500]);
x = [-400:400];
lobf = polyfit(iceprod_est_RCP85_alldata/1e9,YIP*1e3,1);
p = plot(x,lobf(1)*x+lobf(2),'r');
ax2.XTickLabel = []; ylabel('Winter ice production, km^3')

% full scatter plot
ax1 = axes('Position',[0.75,0.15,0.2,0.225]);

YIP = -diceprod_full(:)/1e12;
beta = lm_mrIP_full.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

iceprod_est_full_alldata = beta(1)...
    + beta(2)*dsnow_effect_full(:)...
    + beta(3)*dSIAdeficit_effect_full(:)...
    + beta(4)*dnetdiv_effect_full(:)...
    + beta(5)*dcompdiv_effect_full(:)...
    + beta(6)*dSep_10mT_full(:);

s1 = scatter(iceprod_est_full_alldata/1e9,YIP*1e3); s1.MarkerEdgeColor = c1; s1.Marker = '.'; s1.MarkerEdgeAlpha = 0.4; s1.SizeData = 30;
hold on; grid on; xlim([-500 500]); ylim([-500 500]);
x = [-400:400];
lobf = polyfit(iceprod_est_full_alldata/1e9,YIP*1e3,1);
p = plot(x,lobf(1)*x+lobf(2),'r');

xlabel('linear model solutions, km^3')

txt = text(ax3,10,-400,strcat('R^2 = ',{' '},num2str(round(lm_mrIP_20C.Rsquared.Ordinary,2,'significant')))); txt.FontName = 'Helvetica'; txt.FontSize = 11; txt.Color = c5;
txt = text(ax2,10,-400,strcat('R^2 = ',{' '},num2str(round(lm_mrIP_RCP85.Rsquared.Ordinary,2,'significant')))); txt.FontName = 'Helvetica'; txt.FontSize = 11; txt.Color = c5;
txt = text(ax1,10,-400,strcat('R^2 = ',{' '},num2str(round(lm_mrIP_full.Rsquared.Ordinary,2,'significant')))); txt.FontName = 'Helvetica'; txt.FontSize = 11; txt.Color = c5;

txt = text(ax3,-490,400,strcat('N = ',{' '},num2str(lm_mrIP_20C.NumObservations))); txt.FontName = 'Helvetica'; txt.FontSize = 11; txt.Color = [0.5 0.5 0.5];
txt = text(ax2,-490,400,strcat('N = ',{' '},num2str(lm_mrIP_RCP85.NumObservations))); txt.FontName = 'Helvetica'; txt.FontSize = 11; txt.Color = [0.5 0.5 0.5];
txt = text(ax1,-490,400,strcat('N = ',{' '},num2str(lm_mrIP_full.NumObservations))); txt.FontName = 'Helvetica'; txt.FontSize = 11; txt.Color = [0.5 0.5 0.5];

txa = text(ax,0.12,135,'a)'); txa.FontName = 'Helvetica'; txa.FontSize = 13; txa.FontWeight = 'bold';
txb = text(ax3,-700,350,'b)'); txb.FontName = 'Helvetica'; txb.FontSize = 13; txb.FontWeight = 'bold';
txc = text(ax2,-700,350,'c)'); txc.FontName = 'Helvetica'; txc.FontSize = 13; txc.FontWeight = 'bold';
txd = text(ax1,-700,350,'d)'); txd.FontName = 'Helvetica'; txd.FontSize = 13; txd.FontWeight = 'bold';




ax1.FontName = 'Helvetica'; ax1.FontSize = 12;
ax2.FontName = 'Helvetica'; ax2.FontSize = 12;
ax3.FontName = 'Helvetica'; ax3.FontSize = 12;

set(gcf, 'Color', 'w');
% print(gcf,'-dpdf','~/Documents/MATLAB/project/data/CESM/figures/IceProd_coefficients_bars_uncertainty_all40_scatters.pdf','-bestfit');

%% Reconstruction of individual ensemble members


beta = lm_mrIP_full.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

M = subplot2axes(5,8,0.04,0.01); M(:,2) = M(:,2)+0.1; M(:,1) = M(:,1) + 0.025;
figure;
for n = 1:40

iceprod_est_full_byfull = -mean(iceprod_full(:,n)) + beta(1)...
    + beta(2).*(snow_effect_full(:,n)-mean(snow_effect_full(:,n)))...
    + beta(3).*(SIAdeficit_effect_full(:,n)-mean(SIAdeficit_effect_full(:,n)))...
    + beta(4)*(netdiv_effect_full(:,n)-mean(netdiv_effect_full(:,n)))...
    + beta(5).*(compdiv_effect_full(:,n)-mean(compdiv_effect_full(:,n)))...
    + beta(6)*(Sep_10mT_full(:,n)-mean(Sep_10mT_full(:,n)));

ax=axes('Position',M(n,:));
p = plot(yrs_full, 1e-12*((-1*iceprod_full(:,n))),'color','k'); p.LineWidth = 1;  hold on;
p = plot(yrs_full, 1e-12*(iceprod_est_full_byfull-1*mean(iceprod_full(:,n))-mean(iceprod_est_full_byfull(:))),'color',c1); p.LineWidth = 1; 
text(1930,1.4,num2str(n),'color',c5)
% ax.XTickLabel = [];
ax.YLim = [1.2 2.6]; ax.XLim = [1920 2080]; grid on;
if n==33||n==34||n==35||n==36||n==37||n==38||n==39||n==40
    continue
else
    ax.XTickLabel = []; 
end
end
suptitle('Reconstruction of ice production (x1000 km^3) in individual ensemble members');
% suplabel('Ice production, 1000 km^3','y');
set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/Ice_Production_recon_all_ens_members40.pdf
%% reconstruction with obs
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
X_NOAA = X; lat_SST = lat; lon_SST = lon; Sep_SST_NOAA = Sep_SST_KL; clearvars X lat lon Sep_SST_KL
load snowLG_KL.mat
load KL_Ice_Factories_HadISST.mat X Sep_SST_KL
X_HadISST = X; Sep_SST_HadISST = Sep_SST_KL; clearvars X Sep_SST_KL

% create timeseries

% get the right years!!
% Sep goes with the year following
% Sep_SST and Sep_SIA need to be treated as yrs+1
X_NOAA1 = X_NOAA + 1;
X_HadISST1 = X_HadISST + 1;
X_SIC1 = X_SIC + 1;
% record will be 1983 (first SST year) to 2019 (last snow year)
X = 1983:2019;
fNOAA = find(X_NOAA1==X(1)):find(X_NOAA1==X(end));
fHadISST = find(X_HadISST1==X(1)):find(X_HadISST1==X(end));
f2mT = find(X_2mT==X(1)):find(X_2mT==X(end));
fSIC = find(X_SIC1==X(1)):find(X_SIC1==X(end));
fsnow = find(snowKLM2.yrs==X(1)):find(snowKLM2.yrs==X(end));
fPP = find(X_PP==X(1)):find(X_PP==X(end));

deltaT_obs = -(winter_2mT(f2mT)-273.15+1.8);

total_area = 1.2543e12;
SIAdeficit_effect_obs = (total_area - SIA_sum(fSIC)).*deltaT_obs;
snow_effect_obs = deltaT_obs./snowKLM2.Mean(fsnow)';
compdiv_area = posdiv_area - netdiv_area;
netdiv_effect_obs = deltaT_obs.*netdiv_area(fPP);
posdiv_effect_obs = deltaT_obs.*posdiv_area(fPP);
compdiv_effect_obs = deltaT_obs.*compdiv_area(fPP);

% do reconstruction
XSIAdeficit_effect = dSIAdeficit_effect_full(:)/1e12;
Xnetdiv_effect = dnetdiv_effect_full(:)/1e12;
Xcompdiv_effect = dcompdiv_effect_full(:)/1e12;
Xsnow_effect = dsnow_effect_full(:);
XSep_10mT = dSep_10mT_full(:);

XIP = [Xsnow_effect,XSIAdeficit_effect,Xnetdiv_effect,Xcompdiv_effect,XSep_10mT]; 
YIP = -diceprod_full(:)/1e12;

lm_mrIP_full = fitlm(XIP,YIP); 
beta = lm_mrIP_full.Coefficients.Estimate; beta(3:5) = beta(3:5)/1e12; beta = beta*1e12;

%With NOAA SST
iceprod_est_obs_v1 = beta(1)...
    + beta(2)*(snow_effect_obs-mean(snow_effect_obs))...
    + beta(3)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs))...
    + beta(4)*(netdiv_effect_obs-mean(netdiv_effect_obs))...
    + beta(5)*(compdiv_effect_obs-mean(compdiv_effect_obs))...
    + beta(6)*(Sep_SST_NOAA(fNOAA)-mean(Sep_SST_NOAA(fNOAA)));
%With HadISST
iceprod_est_obs_v2 = beta(1)...
    + beta(2)*(snow_effect_obs-mean(snow_effect_obs))...
    + beta(3)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs))...
    + beta(4)*(netdiv_effect_obs-mean(netdiv_effect_obs))...
    + beta(5)*(compdiv_effect_obs-mean(compdiv_effect_obs))...
    + beta(6)*(Sep_SST_HadISST(fHadISST)-mean(Sep_SST_HadISST(fHadISST)));
figure;
p = plot(X,iceprod_est_obs_v1); p.LineWidth = 2; hold on;
p = plot(X,iceprod_est_obs_v2); p.LineWidth = 2;
%%
M = subplot2axes(1,2,0.075,0.075); M(:,1) = M(:,1) + 0.05;
figure;
ax1 = axes('Position',M(1,:));
p = plot(X, 1e-12*(iceprod_est_obs_v1-iceprod_est_obs_v1(1)),'color',c1); p.LineWidth = 2; hold on; p.LineStyle = ':';
p = plot(X, 1e-12*(iceprod_est_obs_v2-iceprod_est_obs_v2(1)),'color',c1); p.LineWidth = 2; hold on;
% p = plot(X,1e-12*(movmean(iceprod_est_obs_v1,10)-iceprod_est_obs_v1(1))); p.LineWidth = 2; p.Color = c5;
p = plot(X, 1e-12*(beta(2)*(snow_effect_obs-mean(snow_effect_obs)) - beta(2)*(snow_effect_obs(1)-mean(snow_effect_obs)))); p.Color = c2;
p = plot(X, 1e-12*(beta(3)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs)) - beta(3)*(SIAdeficit_effect_obs(1)-mean(SIAdeficit_effect_obs)))); p.Color = c3;
p = plot(X, 1e-12*(beta(4)*(netdiv_effect_obs-mean(netdiv_effect_obs)) - beta(4)*(netdiv_effect_obs(1)-mean(netdiv_effect_obs)))); p.Color = c4;
p = plot(X, 1e-12*(beta(5)*(compdiv_effect_obs-mean(compdiv_effect_obs)) - beta(5)*(compdiv_effect_obs(1)-mean(compdiv_effect_obs)))); p.Color = c7;
p = plot(X, 1e-12*( beta(6)*(Sep_SST_NOAA(fNOAA)-mean(Sep_SST_NOAA(fNOAA))) - beta(6)*(Sep_SST_NOAA(fNOAA(1))-mean(Sep_SST_NOAA(fNOAA))))); p.Color = c6; p.LineStyle = ':';
p = plot(X, 1e-12*( beta(6)*(Sep_SST_HadISST(fHadISST)-mean(Sep_SST_HadISST(fHadISST))) - beta(6)*(Sep_SST_HadISST(fHadISST(1))-mean(Sep_SST_HadISST(fHadISST))))); p.Color = c6;

ylim(ax1,[-0.2 0.3]);
t = title('Full reconstruction'); t.FontWeight = 'normal';
grid on;
% L = legend('Ice production recon.','(1/Snow depth) * \DeltaT','SIA deficit * \DeltaT','Net div * \DeltaT','Comp div * \DeltaT','Sep SST'); L.Location = 'NorthWest'; L.AutoUpdate = 'off';
%%
M(2,1) = M(2,1) - 0.05;
ax2 = axes('Position',M(2,:)); 
p = plot(X,1e-12*(movmean(iceprod_est_obs_v2,10)-iceprod_est_obs_v2(1))); p.LineWidth = 2; p.Color = c1; hold on;
p = plot(X, 1e-12*(beta(2)*(movmean(snow_effect_obs,10)-mean(snow_effect_obs)) - beta(2)*(snow_effect_obs(1)-mean(snow_effect_obs)))); p.LineWidth = 2; p.Color = c2;
p = plot(X, 1e-12*(beta(3)*(movmean(SIAdeficit_effect_obs,10)-mean(SIAdeficit_effect_obs)) - beta(3)*(SIAdeficit_effect_obs(1)-mean(SIAdeficit_effect_obs)))); p.LineWidth = 2; p.Color = c3;
p = plot(X, 1e-12*(beta(4)*(movmean(netdiv_effect_obs,10)-mean(netdiv_effect_obs)) - beta(4)*(netdiv_effect_obs(1)-mean(netdiv_effect_obs)))); p.LineWidth = 2; p.Color = c4;
p = plot(X, 1e-12*(beta(5)*(movmean(compdiv_effect_obs,10)-mean(compdiv_effect_obs)) - beta(5)*(compdiv_effect_obs(1)-mean(compdiv_effect_obs)))); p.LineWidth = 2; p.Color = c7;
p = plot(X, 1e-12*( beta(6)*(movmean(Sep_SST_HadISST(fHadISST),10)-mean(Sep_SST_HadISST(fHadISST))) - beta(6)*(Sep_SST_HadISST(fHadISST(1))-mean(Sep_SST_HadISST(fHadISST))))); p.LineWidth = 2; p.Color = c6;
p = plot([X(1),X(end)],[0 0]); p.LineStyle = ':'; p.Color = 'k';
p = plot(X,1e-12*(movmean(iceprod_est_obs_v1,10)-iceprod_est_obs_v1(1))); p.LineWidth = 2; p.Color = c1; hold on; p.LineStyle = ':';
p = plot(X, 1e-12*(beta(6)*(movmean(Sep_SST_NOAA(fNOAA),10)-mean(Sep_SST_NOAA(fNOAA))) - beta(6)*(Sep_SST_NOAA(fNOAA(1))-mean(Sep_SST_NOAA(fNOAA))))); p.LineWidth = 2; p.Color = c6; p.LineStyle = ':';

ylim([-0.2 0.3]); grid on;
t = suptitle('Reconstructing historical changes in Kara-Laptev ice production'); t.FontName = 'Helvetica'; t.FontWeight = 'bold';
ylabel(ax1,'Ice production, 1000 km^3'); ax2.YTickLabel = []; ax2.XTickLabel = {'','1990','2000','2010','2020'};
t = title('10 year moving means'); t.FontWeight = 'normal';
ax1.FontSize = 13; ax1.FontName = 'Helvetica';
ax2.FontSize = 13; ax2.FontName = 'Helvetica';
set(gcf, 'Color', 'w');
L = legend('Ice production recon.','(1/Snow depth) * \DeltaT','Sep open water * \DeltaT','Net div * \DeltaT','Comp div * \DeltaT','Sep SST'); L.Location = 'best'; L.AutoUpdate = 'off';


%% Reconstructing obs with uncertainty bounds from all ensemble members 

% combine all estimates of regression coefficients
beta_ens_IP_all3 = [beta_ens_IP_20C,beta_ens_IP_RCP85,beta_ens_IP_full];
% this gives us 120 estimates

for i = 1:120
    %With NOAA SST
iceprod_est_obs40_v1(:,i) = beta_ens_IP_all3(1,i)...
    + beta_ens_IP_all3(2,i)*(snow_effect_obs-mean(snow_effect_obs))...
    + beta_ens_IP_all3(3,i)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs))...
    + beta_ens_IP_all3(4,i)*(netdiv_effect_obs-mean(netdiv_effect_obs))...
    + beta_ens_IP_all3(5,i)*(compdiv_effect_obs-mean(compdiv_effect_obs))...
    + beta_ens_IP_all3(6,i)*(Sep_SST_NOAA(fNOAA)-mean(Sep_SST_NOAA(fNOAA)));
%With HadISST
iceprod_est_obs40_v2(:,i) = beta_ens_IP_all3(1,i)...
    + beta_ens_IP_all3(2,i)*(snow_effect_obs-mean(snow_effect_obs))...
    + beta_ens_IP_all3(3,i)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs))...
    + beta_ens_IP_all3(4,i)*(netdiv_effect_obs-mean(netdiv_effect_obs))...
    + beta_ens_IP_all3(5,i)*(compdiv_effect_obs-mean(compdiv_effect_obs))...
    + beta_ens_IP_all3(6,i)*(Sep_SST_HadISST(fHadISST)-mean(Sep_SST_HadISST(fHadISST)));

snow_est_obs40(:,i) = (beta_ens_IP_all3(2,i)*(snow_effect_obs-mean(snow_effect_obs)) - beta_ens_IP_all3(2,i)*(snow_effect_obs(1)-mean(snow_effect_obs)));
SIA_est_obs40(:,i) = (beta_ens_IP_all3(3,i)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs)) - beta_ens_IP_all3(3,i)*(SIAdeficit_effect_obs(1)-mean(SIAdeficit_effect_obs)));
netdiv_est_obs40(:,i) = (beta_ens_IP_all3(4,i)*(netdiv_effect_obs-mean(netdiv_effect_obs)) - beta_ens_IP_all3(4,i)*(netdiv_effect_obs(1)-mean(netdiv_effect_obs))); 
compdiv_est_obs40(:,i) = (beta_ens_IP_all3(5,i)*(compdiv_effect_obs-mean(compdiv_effect_obs)) - beta_ens_IP_all3(5,i)*(compdiv_effect_obs(1)-mean(compdiv_effect_obs)));
SepSST_est_obs40NOAA(:,i) = (beta_ens_IP_all3(6,i)*(Sep_SST_NOAA(fNOAA)-mean(Sep_SST_NOAA(fNOAA))) - beta_ens_IP_all3(6,i)*(Sep_SST_NOAA(fNOAA(1))-mean(Sep_SST_NOAA(fNOAA)))); 
SepSST_est_obs40HadISST(:,i) = (beta_ens_IP_all3(6,i)*(Sep_SST_HadISST(fHadISST)-mean(Sep_SST_HadISST(fHadISST))) - beta_ens_IP_all3(6,i)*(Sep_SST_HadISST(fHadISST(1))-mean(Sep_SST_HadISST(fHadISST))));
end

% take mean and std
msnow_est_obs40 = mean(snow_est_obs40,2); sdsnow_est_obs40 = std(snow_est_obs40,[],2);
mSIA_est_obs40 = mean(SIA_est_obs40,2); sdSIA_est_obs40 = std(SIA_est_obs40,[],2);
mnetdiv_est_obs40 = mean(netdiv_est_obs40,2); sdnetdiv_est_obs40 = std(netdiv_est_obs40,[],2);
mcompdiv_est_obs40 = mean(compdiv_est_obs40,2); sdcompdiv_est_obs40 = std(compdiv_est_obs40,[],2);
mSepSST_est_obs40NOAA = mean(SepSST_est_obs40NOAA,2); sdSepSST_est_obs40NOAA = std(SepSST_est_obs40NOAA,[],2);
mSepSST_est_obs40HadISST = mean(SepSST_est_obs40HadISST,2); sdSepSST_est_obs40HadISST = std(SepSST_est_obs40HadISST,[],2);
miceprod_est_obs40_v1 = mean(iceprod_est_obs40_v1,2); sdiceprod_est_obs40_v1 = std(iceprod_est_obs40_v1,[],2);
miceprod_est_obs40_v2 = mean(iceprod_est_obs40_v2,2); sdiceprod_est_obs40_v2 = std(iceprod_est_obs40_v2,[],2);


% plot 
M = subplot2axes(2,1,0.075,0.075); M(:,1) = M(:,1) + 0.05;
figure;
ax1 = axes('Position',M(1,:));
[h1, hp1] = boundedline(X, 1e-12*(miceprod_est_obs40_v1-miceprod_est_obs40_v1(1)),1e-12*sdiceprod_est_obs40_v1); h1.Color = c1; hp1.FaceColor = c1; hp1.FaceAlpha = 0.5; h1.LineStyle = ':';
[h1, hp1] = boundedline(X, 1e-12*(miceprod_est_obs40_v2-miceprod_est_obs40_v2(1)),1e-12*sdiceprod_est_obs40_v2); h1.Color = c1; hp1.FaceColor = c1; hp1.FaceAlpha = 0.5; h1.LineStyle = '-';
[h1, hp1] = boundedline(X, 1e-12*(msnow_est_obs40-msnow_est_obs40(1)),1e-12*sdsnow_est_obs40); h1.Color = c2; hp1.FaceColor = c2; hp1.FaceAlpha = 0.5; h1.LineStyle = '-';
[h1, hp1] = boundedline(X, 1e-12*(mSIA_est_obs40-mSIA_est_obs40(1)),1e-12*sdSIA_est_obs40); h1.Color = c3; hp1.FaceColor = c3; hp1.FaceAlpha = 0.5; h1.LineStyle = '-';
[h1, hp1] = boundedline(X, 1e-12*(mnetdiv_est_obs40-mnetdiv_est_obs40(1)),1e-12*sdnetdiv_est_obs40); h1.Color = c4; hp1.FaceColor = c4; hp1.FaceAlpha = 0.5; h1.LineStyle = '-';
[h1, hp1] = boundedline(X, 1e-12*(mcompdiv_est_obs40-mcompdiv_est_obs40(1)),1e-12*sdcompdiv_est_obs40); h1.Color = c7; hp1.FaceColor = c7; hp1.FaceAlpha = 0.5; h1.LineStyle = '-';
[h1, hp1] = boundedline(X, 1e-12*(mSepSST_est_obs40NOAA-mSepSST_est_obs40NOAA(1)),1e-12*sdSepSST_est_obs40NOAA); h1.Color = c6; hp1.FaceColor = c6; hp1.FaceAlpha = 0.5; h1.LineStyle = ':';
[h1, hp1] = boundedline(X, 1e-12*(mSepSST_est_obs40HadISST-mSepSST_est_obs40HadISST(1)),1e-12*sdSepSST_est_obs40HadISST); h1.Color = c6; hp1.FaceColor = c6; hp1.FaceAlpha = 0.5; h1.LineStyle = '-';

% p = plot(X, 1e-12*(beta(2)*(snow_effect_obs-mean(snow_effect_obs)) - beta(2)*(snow_effect_obs(1)-mean(snow_effect_obs)))); p.Color = c2;
% p = plot(X, 1e-12*(beta(3)*(SIAdeficit_effect_obs-mean(SIAdeficit_effect_obs)) - beta(3)*(SIAdeficit_effect_obs(1)-mean(SIAdeficit_effect_obs)))); p.Color = c3;
% p = plot(X, 1e-12*(beta(4)*(netdiv_effect_obs-mean(netdiv_effect_obs)) - beta(4)*(netdiv_effect_obs(1)-mean(netdiv_effect_obs)))); p.Color = c4;
% p = plot(X, 1e-12*(beta(5)*(compdiv_effect_obs-mean(compdiv_effect_obs)) - beta(5)*(compdiv_effect_obs(1)-mean(compdiv_effect_obs)))); p.Color = c7;
% p = plot(X, 1e-12*( beta(6)*(Sep_SST_NOAA(fNOAA)-mean(Sep_SST_NOAA(fNOAA))) - beta(6)*(Sep_SST_NOAA(fNOAA(1))-mean(Sep_SST_NOAA(fNOAA))))); p.Color = c6; p.LineStyle = ':';
% p = plot(X, 1e-12*( beta(6)*(Sep_SST_HadISST(fHadISST)-mean(Sep_SST_HadISST(fHadISST))) - beta(6)*(Sep_SST_HadISST(fHadISST(1))-mean(Sep_SST_HadISST(fHadISST))))); p.Color = c6;


t = title('Full reconstruction'); t.FontWeight = 'normal';
grid on; box on;
% L = legend('Ice production recon.','(1/Snow depth) * \DeltaT','SIA deficit * \DeltaT','Net div * \DeltaT','Comp div * \DeltaT','Sep SST'); L.Location = 'NorthWest'; L.AutoUpdate = 'off';

ax2 = axes('Position',M(2,:)); 
% p = plot(yrs_full,-miceprod_full/1e12 - mean(-miceprod_full(63:98)/1e12)+mean(1e-12*(movmean(iceprod_est_obs_v2,10)-iceprod_est_obs_v2(1)))); p.LineWidth =1; p.Color = [0.5 0.5 0.5]; hold on;
p = plot(X,1e-12*(movmean(iceprod_est_obs_v2,10)-iceprod_est_obs_v2(1))); p.LineWidth = 2; p.Color = c1; hold on;
p = plot(X, 1e-12*(beta(2)*(movmean(snow_effect_obs,10)-mean(snow_effect_obs)) - beta(2)*(snow_effect_obs(1)-mean(snow_effect_obs)))); p.LineWidth = 2; p.Color = c2;
p = plot(X, 1e-12*(beta(3)*(movmean(SIAdeficit_effect_obs,10)-mean(SIAdeficit_effect_obs)) - beta(3)*(SIAdeficit_effect_obs(1)-mean(SIAdeficit_effect_obs)))); p.LineWidth = 2; p.Color = c3;
p = plot(X, 1e-12*(beta(4)*(movmean(netdiv_effect_obs,10)-mean(netdiv_effect_obs)) - beta(4)*(netdiv_effect_obs(1)-mean(netdiv_effect_obs)))); p.LineWidth = 2; p.Color = c4;
p = plot(X, 1e-12*(beta(5)*(movmean(compdiv_effect_obs,10)-mean(compdiv_effect_obs)) - beta(5)*(compdiv_effect_obs(1)-mean(compdiv_effect_obs)))); p.LineWidth = 2; p.Color = c7;
p = plot(X, 1e-12*( beta(6)*(movmean(Sep_SST_HadISST(fHadISST),10)-mean(Sep_SST_HadISST(fHadISST))) - beta(6)*(Sep_SST_HadISST(fHadISST(1))-mean(Sep_SST_HadISST(fHadISST))))); p.LineWidth = 2; p.Color = c6;
p = plot([X(1),X(end)],[0 0]); p.LineStyle = ':'; p.Color = 'k';
p = plot(X,1e-12*(movmean(iceprod_est_obs_v1,10)-iceprod_est_obs_v1(1))); p.LineWidth = 2; p.Color = c1; hold on; p.LineStyle = ':';
p = plot(X, 1e-12*(beta(6)*(movmean(Sep_SST_NOAA(fNOAA),10)-mean(Sep_SST_NOAA(fNOAA))) - beta(6)*(Sep_SST_NOAA(fNOAA(1))-mean(Sep_SST_NOAA(fNOAA))))); p.LineWidth = 2; p.Color = c6; p.LineStyle = ':';

 grid on;
t = suptitle('Application of linear model to historical changes'); t.FontName = 'Helvetica'; t.FontWeight = 'bold'; t.Position = t.Position + [0.04 0 0];
% ylabel('Ice production, 1000 km^3'); ylabel(ax1,'Ice production, 1000 km^3');
t = title('10 yr moving means'); t.FontWeight = 'normal';
ax1.FontSize = 13; ax1.FontName = 'Helvetica';
ax2.FontSize = 13; ax2.FontName = 'Helvetica';
ax1.XLim = [1982 2020]; ax2.XLim = [1982 2020]; ax2.YLim = [-0.1 0.2]; ax1.YLim = [-0.35 0.35]; ax1.XTickLabel = []; 
ta = text(ax1,1980.25,0.28,'a)'); ta.FontSize = 13; ta.FontName = 'Helvetica'; ta.FontWeight = 'bold';
tb = text(ax2,1980.25,0.1625,'b)'); tb.FontSize = 13; tb.FontName = 'Helvetica'; tb.FontWeight = 'bold';

set(gcf, 'Color', 'w');
L = legend('Ice production','\beta_1(1/Snow depth)\DeltaT','\beta_2Sep open water\DeltaT','\beta_3Net div\DeltaT','\beta_4Comp div\DeltaT','\beta_5Sep SST'); 
L.Position = [0.14 0.338 0.275 0.210]; L.AutoUpdate = 'off';

sy = suplabel('Winter ice production, 1000 km^3','y'); sy.FontName = 'Helvetica'; sy.FontSize = 13; sy.Position = sy.Position + [0.025 0 0 0];

set(gcf, 'Color', 'w');
% print(gcf,'-dpdf','~/Documents/MATLAB/project/data/CESM/figures/Ice_Production_recon_obsv6_40.pdf','-bestfit');

%%
% find moving means and sum error in quadrature
figure; 
% not sure I want this bit
for n = 1:37
    L = n-4;
    U = n+4;
    if n-4<1
        L = 1;
    end
    if n+4>37
        U = 37;
    end
    v = 1e-12*(miceprod_est_obs40_v1(L:U)-miceprod_est_obs40_v1(1));
    mmiceprod_est_obs40_v1(n) = mean(v);
        v = 1e-12*(sdiceprod_est_obs40_v1(L:U));
    mmerr_iceprod_est_obs40_v1(n) = sqrt(sum(v.^2));
    v = 1e-12*(miceprod_est_obs40_v2(L:U)-miceprod_est_obs40_v1(1));
    mmiceprod_est_obs40_v2(n) = mean(v);
        v = 1e-12*(sdiceprod_est_obs40_v2(L:U));
    mmerr_iceprod_est_obs40_v2(n) = sqrt(sum(v.^2));
    v = 1e-12*(msnow_est_obs40(L:U)-msnow_est_obs40(1));
    mmsnow_est_obs40(n) = mean(v);
        v = 1e-12*(sdsnow_est_obs40(L:U));
    mmerr_snow_est_obs40(n) = sqrt(sum(v.^2));
    v = 1e-12*(mSIA_est_obs40(L:U)-mSIA_est_obs40(1));
    mmSIA_est_obs40(n) = mean(v);
        v = 1e-12*(sdSIA_est_obs40(L:U));
    mmerr_SIA_est_obs40(n) = sqrt(sum(v.^2));
    v = 1e-12*(mnetdiv_est_obs40(L:U)-mnetdiv_est_obs40(1));
    mmnetdiv_est_obs40(n) = mean(v);
        v = 1e-12*(sdnetdiv_est_obs40(L:U));
    mmerr_netdiv_est_obs40(n) = sqrt(sum(v.^2));    
    v = 1e-12*(mcompdiv_est_obs40(L:U)-mcompdiv_est_obs40(1));
    mmcompdiv_est_obs40(n) = mean(v);
        v = 1e-12*(sdcompdiv_est_obs40(L:U));
    mmerr_compdiv_est_obs40(n) = sqrt(sum(v.^2));
    v = 1e-12*(mSepSST_est_obs40NOAA(L:U)-mSepSST_est_obs40NOAA(1));
    mmSepSST_est_obs40NOAA(n) = mean(v);
        v = 1e-12*(sdSepSST_est_obs40NOAA(L:U));
    mmerr_SepSST_est_obs40NOAA(n) = sqrt(sum(v.^2));
    v = 1e-12*(mSepSST_est_obs40HadISST(L:U)-mSepSST_est_obs40HadISST(1));
    mmSepSST_est_obs40HadISST(n) = mean(v);
        v = 1e-12*(sdSepSST_est_obs40HadISST(L:U));
    mmerr_SepSST_est_obs40HadISST(n) = sqrt(sum(v.^2));
end

    
ax2 = axes('Position',M(1,:)); 
% [h1, hp1] = boundedline(X,mmiceprod_est_obs40_v1,mmerr_iceprod_est_obs40_v1); h1.Color = c1; hp1.FaceColor = c1; hp1.FaceAlpha = 0.3; h1.LineStyle = ':'; h1.LineWidth = 2;
[h1, hp1] = boundedline(X,mmiceprod_est_obs40_v2,mmerr_iceprod_est_obs40_v2); h1.Color = c1; hp1.FaceColor = c1; hp1.FaceAlpha = 0.3; h1.LineStyle = '-'; h1.LineWidth = 2;
[h1, hp1] = boundedline(X,mmsnow_est_obs40,mmerr_snow_est_obs40); h1.Color = c2; hp1.FaceColor = c2; hp1.FaceAlpha = 0.3; h1.LineStyle = '-'; h1.LineWidth = 2;
[h1, hp1] = boundedline(X,mmSIA_est_obs40,mmerr_SIA_est_obs40); h1.Color = c3; hp1.FaceColor = c3; hp1.FaceAlpha = 0.3; h1.LineStyle = '-'; h1.LineWidth = 2;
[h1, hp1] = boundedline(X,mmnetdiv_est_obs40,mmerr_netdiv_est_obs40); h1.Color = c4; hp1.FaceColor = c4; hp1.FaceAlpha = 0.3; h1.LineStyle = '-'; h1.LineWidth = 2;
[h1, hp1] = boundedline(X,mmcompdiv_est_obs40,mmerr_compdiv_est_obs40); h1.Color = c7; hp1.FaceColor = c7; hp1.FaceAlpha = 0.3; h1.LineStyle = '-'; h1.LineWidth = 2;
% [h1, hp1] = boundedline(X,mmSepSST_est_obs40NOAA,mmerr_SepSST_est_obs40NOAA); h1.Color = c6; hp1.FaceColor = c6; hp1.FaceAlpha = 0.3; h1.LineStyle = ':'; h1.LineWidth = 2;
[h1, hp1] = boundedline(X,mmSepSST_est_obs40HadISST,mmerr_SepSST_est_obs40HadISST); h1.Color = c6; hp1.FaceColor = c6; hp1.FaceAlpha = 0.3; h1.LineStyle = '-'; h1.LineWidth = 2;
ylim([-0.2 0.3]); grid on;
t = suptitle('Reconstructing historical changes in Kara-Laptev ice production'); t.FontName = 'Helvetica'; t.FontWeight = 'bold';
ylabel('Ice production, 1000 km^3'); ylabel(ax1,'Ice production, 1000 km^3');
t = title('10 year moving means'); t.FontWeight = 'normal';
ax1.FontSize = 13; ax1.FontName = 'Helvetica';
ax2.FontSize = 13; ax2.FontName = 'Helvetica';
ax1.XLim = [1982 2020]; ax2.XLim = [1982 2020]; ax2.YLim = [-0.1 0.25];
set(gcf, 'Color', 'w');
L = legend('Ice production recon.','(1/Snow depth) * \DeltaT','Sep open water * \DeltaT','Net div * \DeltaT','Comp div * \DeltaT','Sep SST'); L.Location = 'best'; L.AutoUpdate = 'off';


set(gcf, 'Color', 'w');
% print(gcf,'-dpdf','~/Documents/MATLAB/project/data/CESM/figures/Ice_Production_recon_obsv5_40.pdf','-bestfit');

%%
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/Ice_Production_recon_obsv4_edit.pdf

% 
% fn = '~/Documents/MATLAB/project/data/CESM/figures/figsmatlab/Ice_Production_recon_obsv4_edit.fig';
% saveas(gcf,fn);
% 
% set(gcf, 'Color', 'w');
% export_fig ~/Documents/MATLAB/project/data/CESM/figures/Ice_Production_recon_obsv2.png -m3
