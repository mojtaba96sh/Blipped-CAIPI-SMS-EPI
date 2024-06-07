function K_sake = SAKE_singleSlice(R1)
%clear
import LmyGhostCorrection.*
import LmyUtility.*
% load_folder=fullfile('.','data');
%load(fullfile(load_folder,'SbData_128kx128ky_lpc'),'epi_kxkyzc_2shot_lpcCor')
% load 'R1.mat';
rawdata = R1;
rawdata = permute(rawdata, [1 3 2]);
R=zeros(size(rawdata));
for j=1:size(rawdata,3)
    for i=1:2:size(rawdata,2)
        R(:,i,j)=rawdata(:,i,j);end
    for i=2:2:size(rawdata,2)
        R(:,i,j)=rawdata(end:-1:1,i,j);end
end
epi_kxkyzc_raw = R;
[epi_kxkyzc_lpcCor,phasepara]=oneDimLinearCorr_entropy(epi_kxkyzc_raw,1);
% [epi_kxkyzc_lpcCor]=oneDimLinearCorr_parameter(epi_kxkyzc_raw,1,phasepara);
epi_kxkyzc_2shot_lpcCor = zeros(size(epi_kxkyzc_lpcCor,1),size(epi_kxkyzc_lpcCor,2),1,size(epi_kxkyzc_lpcCor,3));
epi_kxkyzc_2shot_lpcCor(:,:,1,:) = epi_kxkyzc_lpcCor;
nSlice= 1;
nCoil=size(epi_kxkyzc_2shot_lpcCor,4);
%%
ncalib = 128;
% threshold_list=[4*ones(1,10),4.5*ones(1,15),5*ones(1,30),4.5*ones(1,15),4*ones(1,10)]; % good for phantom
% threshold_list=[linspace(4,5,40),linspace(5,4,40)];
% threshold_list=[linspace(3,4,5),linspace(4,4,15)];
threshold_list=3;
ksize = [3,3]; % ESPIRiT kernel-window-size
sakeIter = 2;
% wnthresh = 4.5; % 3 or 4 good for brain

epi_kxkyzc_2shot_sakeCor=zeros(ncalib,ncalib,nSlice,nCoil);
epi_kxkyzc_2shot_sakeCor_fullCoils=cat(4,epi_kxkyzc_2shot_sakeCor,epi_kxkyzc_2shot_sakeCor);
for iSlice=1:1
    % for iSlice=5
    disp(iSlice)
    DATA=squeeze(epi_kxkyzc_2shot_lpcCor(:,:,iSlice,:));
    nSeg=2;
    % convert shot to VCC
    DATA_org=DATA;
    for iSeg=2:nSeg
        DATA=cat(3,DATA,circshift(DATA_org,-iSeg+1,2));
        % DATA=cat(3,DATA,circshift(DATA_org,0,2));
    end
    if exist('threshold_list','var')
        wnthresh = threshold_list(iSlice); % Window-normalized number of singular values to threshold
    end
    [sx,sy,Nc] = size(DATA);
    mask = zeros(size(DATA,1),size(DATA,2));
    mask(:,1:nSeg:end) = 1;
    % mask = ones(size(mask));

    DATA2recon=DATA.* repmat(mask,[1,1,size(DATA,3)]);
    DATAc = DATA;
    calibc = crop_new(DATAc,[ncalib,ncalib,size(DATA,3)]);
    
    %% Perform SAKE reconstruction to recover the calibration area
    im = ifft2c(DATAc);
    disp('Performing SAKE recovery of calibration');
    tic; calib_sake = SAKEwithInitialValue(calibc, [ksize], wnthresh,sakeIter, 0,repmat(crop_new(mask,[ncalib,ncalib]),[1 1 size(DATA,3)]));toc
    % calib_sake(:,:,end/4+1:end)=circshift(calib_sake(:,:,end/4+1:end),[0 1])

    calib_sake(:,:,end/2+1:end)=circshift(calib_sake(:,:,end/2+1:end),[0 1]);
    % calib_sake(:,:,end/2+1:end)=circshift(calib_sake(:,:,end/2+1:end),[0 0]);

    % calib_sake(:,:,end*3/4+1:end)=circshift(calib_sake(:,:,end*3/4+1:end),[0 1]);
    epi_kxkyzc_2shot_sakeCor_fullCoils(:,:,iSlice,:)= calib_sake;
    % a=ifft2c(calib_sake(:,:,1:end/4))+ifft2c(calib_sake(:,:,end/4+1:end/2));
    a=ifft2c(calib_sake(:,:,1:end/2))+ifft2c(calib_sake(:,:,end/2+1:end));
    % b=ifft2c(calib_sake(:,:,end/2+1:end*3/4))+ifft2c(calib_sake(:,:,end*3/4+1:end));
    % epi_kxkyzc_2shot_sakeCor(:,:,iSlice,:)=fft2c(a+b)/4;
    epi_kxkyzc_2shot_sakeCor(:,:,iSlice,:)=fft2c(a)/2;
end
disp('Done')
epi_kxkyzc_2shot_sakeCor_new = squeeze(epi_kxkyzc_2shot_sakeCor(:,:,1,:));
image_sake = (sum(abs(ifft2c(squeeze(epi_kxkyzc_2shot_sakeCor))).^2,3)).^0.5; image_sake = image_sake/max(image_sake(:));

image_lpc = (sum(abs(ifft2c(squeeze(epi_kxkyzc_2shot_lpcCor))).^2,3)).^0.5; image_lpc = image_lpc/max(image_lpc(:));
% figure(); subplot(121); imshow(rot90(image_sake,1),[0 1]); title('SAKE'); subplot(122); imshow(rot90(image_lpc,1),[0 1]); title('LPC')

figure(); imshow(rot90(image_sake,1),[0 1]); title('SAKE'); 

K_sake = epi_kxkyzc_2shot_sakeCor_new; 


end