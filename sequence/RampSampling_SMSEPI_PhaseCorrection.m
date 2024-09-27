
% phase correction and reconstruction of ramp sampled SMS epi data using navigator data and my linear phase correction function

%% load data
% path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\mb2_brain_qingping_aug24';
% path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\high slew rate sep2024\brain mb4';
% path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\high slew rate sep2024\Ariel_last_mb2\CimaX_20ch';
% % path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\high slew rate sep2024\Ariel_last_mb2\CimaX_64ch';
% path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\high slew rate sep2024\Ariel_last_mb4\CimaX_20ch';
% path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\structural phantom\mb4';
% path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\spherical phantom\mb2';
% path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\acr phantom(second)\mb4';
% path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\siemens spherical tubes\mb2';
% path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\spherical phantom\mb2';
% path = 'E:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\siemens spherical tubes\mb2';
% path = 'E:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\siemens spherical tubes\mb2_notMain';
% path = 'E:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\siemens spherical tubes\last_mb4';
% path = 'E:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\siemens spherical tubes\lassttt_mb4';
path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\siemens spherical tubes\lassttt_mb4';
path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\high slew rate sep2024\Ariel_last_mb2\CimaX_20ch';
path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\high slew rate sep2024\Ariel_last_mb4\CimaX_20ch';
path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\spherical phantom\mb2';
path = 'D:\MojiProject\caipismsepi_data4esmrmb_presentation_phantoms\spherical phantom\mb4';
path = 'G:\MojiProject\Prisma Measurments\spherical phantom\coil 64\mb4';
path = 'G:\MojiProject\Prisma Measurments\structural phantom\coil 64\mb4';
path = 'D:\MojiProject\Prisma Measurments\spherical phantom\coil 64\mb4';
path = 'D:\MojiProject\Prisma Measurments\brain\coil 64\mb4';
path = 'D:\MojiProject\Prisma Measurments\spherical phantom\coil 64\mb4';
path = 'D:\MojiProject\Prisma Measurments\brain\coil 20\mb4';
path = 'D:\MojiProject\Prisma Measurments\brain\coil 64\mb4';
path = 'C:\Users\ladm-shafiekh\Desktop\ESMRMB Presentation and Results\Data and Seq\Structural Phantom\mb2';
path = 'C:\Users\ladm-shafiekh\Desktop\ESMRMB Presentation and Results\Data and Seq\Structural Phantom\mb4';
path = 'C:\Users\ladm-shafiekh\Desktop\ESMRMB Presentation and Results\Data and Seq\Spherical Phantom\mb2';
path = 'C:\Users\ladm-shafiekh\Desktop\ESMRMB Presentation and Results\Data and Seq\Spherical Phantom\mb4';
path = 'C:\Users\ladm-shafiekh\Desktop\ESMRMB Presentation and Results\Data and Seq\Brain\mb2';
path = 'C:\Users\ladm-shafiekh\Desktop\ESMRMB Presentation and Results\Data and Seq\Brain\mb4';
% path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\optimized_smsepi\mb2';
nF = 1;
[rawdata_epi, ktraj_adc_epi, ~, ~, ~] = Read_RampSampled_EPI_Data(nF, path);

%%
nF = 2;
[rawdata_calib, ktraj_adc_calib, ~, ~, ~] = Read_RampSampled_EPI_Data(nF, path);


%% first use the epi_ramp2cart.m (using interp1 matlab function) to change the ramp sampled data to cartesian data

Ny = 128;
[cartData_sms] = epi_ramp2cart(rawdata_epi,ktraj_adc_epi,Ny);
[calData_sms] = epi_ramp2cart(rawdata_calib,ktraj_adc_calib,Ny);

% [cartData_sms,~] = epi_ramp2cart_nufft(rawdata_epi,ktraj_adc_epi,Ny);
% [calData_sms,~] = epi_ramp2cart_nufft(rawdata_calib,ktraj_adc_calib,Ny);

% [cartData_sms,xky_epi] = epi_ramp2cart_nufft(rawdata_epi,ktraj_adc_epi,Ny);
% [calData_sms,xky_calib] = epi_ramp2cart_nufft(rawdata_calib,ktraj_adc_calib,Ny);

% [cartData_sms] = epi_ramp2cart_version2(rawdata_epi,ktraj_adc_epi,Ny);
% [calData_sms] = epi_ramp2cart_version2(rawdata_calib,ktraj_adc_calib,Ny);

%% second use the Linear_Phase_Correction_Concatenate_withoutInvert.m to correct the phase

[PhaseCorrected_rawdata,Image] = Linear_Phase_Correction_Concatenate_withoutInvert(cartData_sms,calData_sms);

% [PhaseCorrected_rawdata,Image] = Linear_Phase_Correction_xky_nufftF_zappas(xky_epi,xky_calib);

% [PhaseCorrected_rawdata,Image] = Linear_Phase_Correction_xky_nufftF(xky_epi,xky_calib);