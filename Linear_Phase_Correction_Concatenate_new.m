function [PhaseCorrected_rawdata,Image] = Linear_Phase_Correction_Concatenate_new(rawdata,calibdata)
c = size(rawdata,2);
Nphase = size(rawdata,3);
Nread = size(rawdata,1);
%PhaseCorrected_rawdata = zeros(c*Nread,Nphase);
%PhaseCorrected_rawdata_hybrid = zeros(c*Nread,Nphase);
%Image = zeros(size(rawdata));
%CorrectionMat = zeros(size(rawdata));
R = zeros(Nread,Nphase,c);
RC = zeros(Nread,Nphase,c);
RC_invert=zeros(size(RC));
R_invert=zeros(size(R));
for j=1:c
    R(:,:,j)=rawdata(:,j,:);
    RC(:,:,j)=calibdata(:,j,:);
end
for j=1:c
    for i=1:2:Nphase
        RC_invert(:,i,j)=RC(:,i,j);end
    for i=2:2:Nphase
        RC_invert(:,i,j)=RC(end:-1:1,i,j);end
    for i=1:2:Nphase
        R_invert(:,i,j)=R(:,i,j);end
    for i=2:2:Nphase
        R_invert(:,i,j)=R(end:-1:1,i,j);end
end
i1RC_invert = zeros(Nread,Nphase,c);
i1R_invert = zeros(Nread,Nphase,c);
for j=1:c
    i1RC_invert(:,:,j) = ifftshift(ifft(ifftshift(squeeze(RC_invert(:,:,j)),1),[],1),1);
    i1R_invert(:,:,j) = ifftshift(ifft(ifftshift(squeeze(R_invert(:,:,j)),1),[],1),1);    
end
i1RC_invert_new = zeros(c*Nread,Nphase);
% i1R_invert_new = zeros(c*Nread,Nphase);
for j=1:c
%     i1R_invert_new((j-1)*Nread+1:j*Nread,:) = squeeze(i1R_invert(:,:,j));
    i1RC_invert_new((j-1)*Nread+1:j*Nread,:) = squeeze(i1RC_invert(:,:,j));
end

%R_invert_new = zeros(c*Nread,Nphase);
%RC_invert_new = zeros(c*Nread,Nphase);
%for j=1:c
%    R_invert_new((j-1)*Nread+1:j*Nread,:) = squeeze(R_invert(:,:,j));
%    RC_invert_new((j-1)*Nread+1:j*Nread,:) = squeeze(RC_invert(:,:,j));
%end
%i1RC_invert = ifftshift(ifft(ifftshift(RC_invert_new,1),[],1),1);
%i1R_invert = ifftshift(ifft(ifftshift(R_invert_new,1),[],1),1);

C=zeros(c*Nread,Nphase-1);
for i=1:Nphase-1
    if mod(i,2)==1
        C(:,i)=i1RC_invert_new(:,i+1)./i1RC_invert_new(:,i);
    else
        C(:,i)=i1RC_invert_new(:,i)./i1RC_invert_new(:,i+1);
    end
end
CU=zeros(c*Nread,Nphase-1);
for j=1:c
    CU((j-1)*Nread+Nread/2:-1:1+(j-1)*Nread,:) = unwrap(angle(C((j-1)*Nread+Nread/2:-1:1+(j-1)*Nread,:)));
    CU((j-1)*Nread+Nread/2:Nread+(j-1)*Nread,:) = unwrap(angle(C((j-1)*Nread+Nread/2:Nread+(j-1)*Nread,:)));
end
% for j=1:c
%     CU((j-1)*Nread+Nread/2:-1:1+(j-1)*Nread,:) = abs(C((j-1)*Nread+Nread/2:-1:1+(j-1)*Nread,:));
%     CU((j-1)*Nread+Nread/2:Nread+(j-1)*Nread,:) = abs(C((j-1)*Nread+Nread/2:Nread+(j-1)*Nread,:));
% end
thereshold = 0.3;
% mask = Threshold(i1RC_invert_new,c,thereshold);
degree = 2;
mask = new_Th(CU,c,thereshold,degree);
% thereshold = 0.1;
% thereshold_2 = 0.4;
% mask = new_Threshold(i1RC_invert_new,c,thereshold,thereshold_2);
% L1 = zeros(1,c);
% L2 = zeros(1,c);
% for j=1:c
%     [L1(j),L2(j)] = Threshold_Phase(CU((j-1)*Nread+1:j*Nread,:),thereshold);
% end
% L1 = floor(mean(L1));
% L2 = floor(mean(L2));
% % % % % % % % CU_new = zeros(size(CU));
% % % % % % % % for j=1:c
% % % % % % % %     CU_new((j-1)*Nread+L1:(j-1)*Nread+L2,:)=CU((j-1)*Nread+L1:(j-1)*Nread+L2,:);
% % % % % % % % end
Cest = zeros(size(CU));
% x1 = [zeros(L1-1,1);(1:L2-L1+1)';zeros(Nread-L2,1)];
% x2 = [zeros(L1-1,1);ones(L2-L1+1,1);zeros(Nread-L2,1)];
x = repmat([(1:Nread)',ones(Nread,1)],c,1);
xm=x.*mask;
ab= zeros(2,Nphase-1);
for i=1:Nphase-1
    %Cest(:,i) = x*inv(x'*x)*(x'*CU(:,i));
    ab(:,i)=inv(xm'*xm)*(xm'*CU(:,i));
end
ab=(ab(:,1)+ab(:,2))/2;
%ab=(ab(:,1)+ab(:,2)+ab(:,3))/3;
Cest1 = x(1:Nread,:)*ab;
%Cest = zeros(c*Nread,Nphase);
%Cest(:,1:Nphase-1) = Cest1;
D=exp(-1i*Cest1);
PhaseCorrected_rawdata_hybrid=i1R_invert;
% for j=1:c
%     D_new = D((j-1)*Nread+1:j*Nread);
%     PhaseCorrected_rawdata_hybrid(:,2:2:end,j)=i1R_invert(:,2:2:end,j).*D_new(:,ones(Nphase/2,1),ones(1,1));
% end

PhaseCorrected_rawdata_hybrid(:,2:2:end,:)=i1R_invert(:,2:2:end,:).*D(:,ones(Nphase/2,1),ones(c,1));
%PhaseCorrected_rawdata = zeros(size(rawdata));
PhaseCorrected_rawdata=fftshift(fft(fftshift(PhaseCorrected_rawdata_hybrid,1),[],1),1);
Kzp = zeros(3*Nread,3*Nphase,c);
Kzp(Nread+1:2*Nread,Nphase+1:2*Nphase,:) = PhaseCorrected_rawdata;
% Image = fftshift(ifft(fftshift(PhaseCorrected_rawdata_hybrid,2),[],2),2);
Image = zeros(size(Kzp));
for i=1:c
    Image(:,:,i) = fftshift(ifft2(ifftshift(squeeze(Kzp(:,:,i)))));
end
Image(find(isnan(Image))) = 0;
ImageSoS=sum(abs(Image).^2,3).^0.5;
figure(); imshow(rot90(squeeze(ImageSoS),1),[]);
%figure();imshow(abs(squeeze(Image(:,j,:))),[]);
%end
% Image_channels = zeros(Nread,c,Nphase);
% for j=1:c
%     Image_channels(:,j,:) = Image((j-1)*Nread+1:j*Nread,:);
%     PhaseCorrected_rawdata(:,j,:) = PhaseCorrected_rawdata_hybrid((j-1)*Nread+1:j*Nread,:);
%     PhaseCorrected_rawdata(:,j,:) = fftshift(fft(fftshift(squeeze(PhaseCorrected_rawdata(:,j,:)),1),[],1),1);
% end
% FinalImage=sqrt(squeeze(sum((abs(Image_channels)).^2,2)));
% figure();imshow(FinalImage,[]);
end






%%
% brain data niklas
% S1
% mask(611:636,:)=1;mask(746:777,:)=1;mask(882:908,:)=1;mask(1152:1179,:)=1;mask(1292:1326,:)=1;mask(1561:1594,:)=1;mask(2512:2542,:)=1;mask(2648:2668,:)=1;mask(3329:3355,:)=1;mask(3463:3498,:)=1;mask(3584:3610,:)=1;mask(5226:5242,:)=1;mask(5499:5528,:)=1;mask(5640:5657,:)=1;mask(6046:6064,:)=1;mask(6183:6202,:)=1;mask(6314:6335,:)=1;

