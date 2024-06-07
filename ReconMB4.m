

% load 'R.mat'; load 'Rc_GzPrewind.mat'; load 'R1.mat'; load 'Rc1.mat'; load 'R2.mat'; load 'Rc2.mat'; 
% load 'R3.mat'; load 'Rc3.mat'; load 'R4.mat'; load 'Rc4.mat'; 
mb = 4;
% R = permute(R,[1 3 2]);
% Rc = permute(Rc,[1 3 2]);
% R1 = permute(R1,[1 3 2]);
% Rc1 = permute(Rc1,[1 3 2]);
% R2 = permute(R2,[1 3 2]);
% Rc2 = permute(Rc2,[1 3 2]);
[KK,Image] = Linear_Phase_Correction_Concatenate_new(R,Rc_1);
[K1,Image1] = Linear_Phase_Correction_Concatenate_new(R1,Rc1_1);
[K2,Image2] = Linear_Phase_Correction_Concatenate_new(R2,Rc2_1);
[K3,Image3] = Linear_Phase_Correction_Concatenate_new(R3,Rc3_1);
[K4,Image4] = Linear_Phase_Correction_Concatenate_new(R4,Rc4_1);
C1 = zeros(size(KK));
C2 = C1;
C3 = C1;
C4 = C1;

for j=-(size(C1,2))/2:size(C1,2)/2-1
    C1(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*(0));
    C2(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*(1));
    C3(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*(2));
    C4(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*(3));
end
slice_Grappa4



Phi = zeros(size(K1));
fov = [256, 256, 180]*1e-3;
Nz = 64;
mb = 4;
sliceThickness = fov(3)/Nz;            % slice
slicesPerBand=Nz/mb;
bandStep = sliceThickness*slicesPerBand/2;
Ablip = 1/bandStep/mb;
Zoffset = bandStep;
phi = mb*pi*Ablip*Zoffset;
Nphase = size(K1,2);
for a=-Nphase/2:2:Nphase/2-1
    % Phi(:,a+Nphase/2+1,:) = exp(-1i*phi/2*a);
    Phi(:,a+Nphase/2+1,:) = exp(-1i*phi/2);
end
for a=-Nphase/2+1:2:Nphase/2
    % Phi(:,a+Nphase/2+1,:) = exp(1i*phi/2*a);
    Phi(:,a+Nphase/2+1,:) = exp(1i*phi/2);
end

KK = KK.*Phi;
K1 = K1.*Phi;
K2 = K2.*Phi;