% this function gets the ramp sampled EPI and navigator data and then uses 1d nufft for each readout line separately and 
% then we use the linear phase correction (first we need to adapt it with this new input)

%% get EPI and navigator data
path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\mb2';
% EPI data
nF = 6;
[EPI_data, ktraj_adc_EPI, ADClen_EPI, Nread_EPI, Nch_EPI] = Read_RampSampled_EPI_Data(nF, path);
% navigator data
nF = 7;
[calib_data, ktraj_adc_calib, ADClen_calib, Nread_calib, Nch_calib] = Read_RampSampled_EPI_Data(nF, path);

%% apply 1d nufft in readout direction
Nread = 128;
Nphase = Nread_EPI;
xky_EPI = zeros(Nread,Nphase,Nch_EPI);
xky_calib = zeros(Nread,Nphase,Nch_EPI);
for a=1:Nphase
    for c=1:Nch_EPI
        % DCF for EPI data
        Kxx_epi = zeros(1,ADClen_EPI+2);
        Kxx_epi(1) = ktraj_adc_EPI(1,(a-1)*ADClen_EPI+1)-(ktraj_adc_EPI(1,(a-1)*ADClen_EPI+2)-ktraj_adc_EPI(1,(a-1)*ADClen_EPI+1));
        Kxx_epi(end) = ktraj_adc_EPI(1,a*ADClen_EPI)+(ktraj_adc_EPI(1,a*ADClen_EPI)-ktraj_adc_EPI(1,a*ADClen_EPI-1));
        Kxx_epi(2:end-1) = ktraj_adc_EPI(1,(a-1)*ADClen_EPI+1:a*ADClen_EPI);
        D_epi = zeros(1,length(Kxx_epi)-1);
        for i=1:length(Kxx_epi)-1
            D_epi(i) = (Kxx_epi(i)+Kxx_epi(i+1))/2;
        end
        DCF_EPI = diff(D_epi);

        % DCF for calib data
        Kxx_calib = zeros(1,ADClen_calib+2);
        Kxx_calib(1) = ktraj_adc_calib(1,(a-1)*ADClen_calib+1)-(ktraj_adc_calib(1,(a-1)*ADClen_calib+2)-ktraj_adc_calib(1,(a-1)*ADClen_calib+1));
        Kxx_calib(end) = ktraj_adc_calib(1,a*ADClen_calib)+(ktraj_adc_calib(1,a*ADClen_calib)-ktraj_adc_calib(1,a*ADClen_calib-1));
        Kxx_calib(2:end-1) = ktraj_adc_calib(1,(a-1)*ADClen_calib+1:a*ADClen_calib);
        D_calib = zeros(1,length(Kxx_calib)-1);
        for i=1:length(Kxx_calib)-1
            D_calib(i) = (Kxx_calib(i)+Kxx_calib(i+1))/2;
        end
        DCF_calib = diff(D_calib);

        % apply 1d nufft for each readout line of epi data separately
        xky_EPI(:,a,c) = (FGG_1d_type1(squeeze(EPI_data((a-1)*ADClen_EPI+1:a*ADClen_EPI,1,c)).*DCF_EPI',ktraj_adc_EPI(1,(a-1)*ADClen_EPI+1:a*ADClen_EPI),Nread,12));

        % apply 1d nufft for each readout line of calib data separately
        xky_calib(:,a,c) = (FGG_1d_type1(squeeze(calib_data((a-1)*ADClen_calib+1:a*ADClen_calib,1,c)).*DCF_calib',ktraj_adc_calib(1,(a-1)*ADClen_calib+1:a*ADClen_calib),Nread,12));
    end
end

%% use revised linear phase correction
% [PhaseCorrected_rawdata,Image] = Linear_Phase_Correction_revised_4_rampEPI_1dNUFFT(xky_EPI,fftshift(xky_calib,2));
[PhaseCorrected_rawdata,Image] = Linear_Phase_Correction_revised_4_rampEPI_1dNUFFT(xky_EPI,xky_calib);






