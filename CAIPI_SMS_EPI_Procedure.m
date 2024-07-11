%% load SMS data and its calibration data as well as single epi data of all slices
% load R.mat;
% load Rc.mat;
% load R1.mat;
% load R2.mat;
% load R3.mat;
% load R4.mat;

%% Ghost correction: for SMS data we use LPC and for single EPI data we use Antropy-based LPC and SAKE
[KK,~] = Linear_Phase_Correction_Concatenate_new(R,Rc);
K1 = SAKE_singleSlice(R1);
K2 = SAKE_singleSlice(R2);
K3 = SAKE_singleSlice(R3);
K4 = SAKE_singleSlice(R4);


%% Reconstruction by slice-GRAPPA
% mb = 2; % multi band factor
% Images = slicee_GRAPPA_mb2(KK,K1,K2);

mb = 4; 
Images = slicee_GRAPPA_mb4(KK,K1,K2,K3,K4);
