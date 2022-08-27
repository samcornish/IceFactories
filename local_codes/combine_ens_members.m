% combine variables from ensemble members 1:35 and 36:40 into a common file
clear
close all

addpath ~/Documents/MATLAB/project/data/CESM/

load CESM_iceprod_vars_35.mat
Fw_20C_35 = Fw_20C;
Fw_RCP85_35 = Fw_RCP85;
GR_20C_35 = GR_20C;
GR_RCP85_35 = GR_RCP85;
iceprod_20C_35 = iceprod_20C;
iceprod_RCP85_35 = iceprod_RCP85;
net_areadiv_20C_35 = net_areadiv_20C;
net_areadiv_RCP85_35 = net_areadiv_RCP85;
pos_areadiv_20C_35 = pos_areadiv_20C;
pos_areadiv_RCP85_35 = pos_areadiv_RCP85;
hs_20C_35 = hs_20C;
hs_RCP85_35 = hs_RCP85;
SepSIA_20C_35 = SepSIA_20C;
SepSIA_RCP85_35 = SepSIA_RCP85;
SAT_20C_35 = SAT_20C;
SAT_RCP85_35 = SAT_RCP85;
frz_area_days_20C_35 = frz_area_days_20C;
frz_area_days_RCP85_35 = frz_area_days_RCP85;

load CESM_iceprod_vars_36to40.mat

Fw_20C(:,1:35) = Fw_20C_35;
Fw_RCP85(:,1:35) = Fw_RCP85_35;
GR_20C(:,1:35) = GR_20C_35;
GR_RCP85(:,1:35) = GR_RCP85_35;
iceprod_20C(:,1:35) = iceprod_20C_35;
iceprod_RCP85(:,1:35) = iceprod_RCP85_35;
net_areadiv_20C(:,1:35) = net_areadiv_20C_35;
net_areadiv_RCP85(:,1:35) = net_areadiv_RCP85_35;
pos_areadiv_20C(:,1:35) = pos_areadiv_20C_35;
pos_areadiv_RCP85(:,1:35) = pos_areadiv_RCP85_35;
hs_20C(:,1:35) = hs_20C_35;
hs_RCP85(:,1:35) = hs_RCP85_35;
SepSIA_20C(:,1:35) = SepSIA_20C_35;
SepSIA_RCP85(:,1:35) = SepSIA_RCP85_35;
SAT_20C(:,1:35) = SAT_20C_35;
SAT_RCP85(:,1:35) = SAT_RCP85_35;
frz_area_days_20C(:,1:35) = frz_area_days_20C_35;
frz_area_days_RCP85(:,1:35) = frz_area_days_RCP85_35;

save ~/Documents/MATLAB/project/data/CESM/CESM_iceprod_vars_40.mat GR_20C GR_RCP85 iceprod_20C iceprod_RCP85 net_areadiv_20C net_areadiv_RCP85 pos_areadiv_20C pos_areadiv_RCP85 hs_20C hs_RCP85 SepSIA_20C SepSIA_RCP85 SAT_20C SAT_RCP85 frz_area_days_20C frz_area_days_RCP85 Fw_20C Fw_RCP85

%% monthly I2O

clear
close all

addpath ~/Documents/MATLAB/project/data/CESM/

load CESM_monthly_I2O.mat
monthly_I2O_20C_35 = monthly_I2O_20C;
monthly_iceprod_20C_35 = monthly_iceprod_20C;
monthly_melt_20C_35 = monthly_melt_20C;
monthly_I2O_RCP85_35 = monthly_I2O_RCP85;
monthly_iceprod_RCP85_35 = monthly_iceprod_RCP85;
monthly_melt_RCP85_35 = monthly_melt_RCP85;

load CESM_monthly_I2O_36to40.mat
monthly_I2O_20C(:,:,1:35) = monthly_I2O_20C_35;
monthly_iceprod_20C(:,:,1:35) = monthly_iceprod_20C_35;
monthly_melt_20C(:,:,1:35) = monthly_melt_20C_35;
monthly_I2O_RCP85(:,:,1:35) = monthly_I2O_RCP85_35;
monthly_iceprod_RCP85(:,:,1:35) = monthly_iceprod_RCP85_35;
monthly_melt_RCP85(:,:,1:35) = monthly_melt_RCP85_35;

save ~/Documents/MATLAB/project/data/CESM/CESM_monthly_I2O_40.mat