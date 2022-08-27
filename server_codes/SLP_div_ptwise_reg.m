% script to compute the pointwise regression of mean winter SLP against the mean positive divergence of that winter

clear
close all

addpath /home/ocean_personal_data/samc/CESM/variables/

load CESM_div_I2O_reg_20C_RCP85.mat meanposdiv_20C meanposdiv_RCP85
load CESM_SLP.mat

% no deseasonalising required as we are dealing with winter means

% reshape pressure field into array with one spatial dimension
sz_p = size(mean_slpE_20C);
p_20C = reshape(mean_slpE_20C,sz_p(1)*sz_p(2),sz_p(3),sz_p(4));

sz_p = size(mean_slpE_RCP85);
p_RCP85 = reshape(mean_slpE_RCP85,sz_p(1)*sz_p(2),sz_p(3),sz_p(4));

p_20C1980 = p_20C(:,1:60,:);
p_vec = reshape(p_20C1980,sz_p(1)*sz_p(2),60*sz_p(4));
meanposdiv_1980 = meanposdiv_20C(1:60,:);
mpd_vec = reshape(meanposdiv_1980,60*sz_p(4),1);
% regressions
sp = size(p_vec);
for i = 1:sp(1)
	X = [ones(sp(2),1),p_vec(i,:)'/1e5]; % divide by 1e5 to get numbers on same order of magnitude
	Y = mpd_vec;
	[b,bint,r,rint,stats] = regress(Y,X);
	B0_vec(i) = b(1); B1_vec(i) = b(2); R2_vec(i) = stats(1); P_vec(i) = stats(4);
	fprintf('finished pt. # %d\n',i)
end

% reshape back into spatial array
B0 = reshape(B0_vec,sz_p(1),sz_p(2));
B1 = reshape(B1_vec,sz_p(1),sz_p(2))/1e5;
R2 = reshape(R2_vec,sz_p(1),sz_p(2));
Pval = reshape(P_vec,sz_p(1),sz_p(2));

[val ind] = max(R2(:));

[lowi lowj] = ind2sub([sz_p(1),sz_p(2)],ind);

lowP_20C = squeeze(mean_slpE_20C(lowi,lowj,:,:));
lowP_RCP85 = squeeze(mean_slpE_RCP85(lowi,lowj,:,:));

mlowP_20C = mean(lowP_20C,2);
mlowP_RCP85 = mean(lowP_RCP85,2);

meanposdiv_20C_est = B0(lowi,lowj) + B1(lowi,lowj)*mlowP_20C;
meanposdiv_RCP85_est = B0(lowi,lowj) + B1(lowi,lowj)*mlowP_RCP85;
 



%% multiple regression using Dec SIV as well as Barents Low to predict mean pos div

load CESM_SIT_SIC_SIV_20C_RCP85.mat

Dec_SIV_20C = squeeze(monthly_SIV_20C(:,3,:));
Dec_SIV_RCP85 = squeeze(monthly_SIV_RCP85(:,3,:));
SIV = reshape(Dec_SIV_20C(1:60,:),60*32,1)/1e12;
LOW = reshape(lowP_20C(1:60,:),60*32,1)/1e5;
X = [ones(length(SIV),1),SIV,LOW]; Y = mpd_vec;
[b,bint,r,rint,stats] = regress(Y,X);
B0_mr = b(1); B1_mr = b(2)/1e12; B2_mr = b(3)/1e5; R2_mr = stats(1);

save /home/ocean_personal_data/samc/CESM/variables/slp_mpd_regression.mat B0 B1 R2 Pval lat lon mlowP_20C mlowP_RCP85 lowP_20C lowP_RCP85 meanposdiv_20C_est meanposdiv_RCP85_est B0_mr B1_mr B2_mr R2_mr

