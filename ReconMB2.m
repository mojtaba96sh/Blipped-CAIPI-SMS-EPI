

load 'R.mat'; load 'Rc.mat'; load 'R1.mat'; load 'Rc1.mat'; load 'R2.mat'; load 'Rc2.mat'; 
mb = 2;
% R = permute(R,[1 3 2]);
% Rc = permute(Rc,[1 3 2]);
% R1 = permute(R1,[1 3 2]);
% Rc1 = permute(Rc1,[1 3 2]);
% R2 = permute(R2,[1 3 2]);
% Rc2 = permute(Rc2,[1 3 2]);
[KK,Image] = Linear_Phase_Correction_Concatenate_new(R,Rc);
[K1,Image1] = Linear_Phase_Correction_Concatenate_new(R1,Rc1);
[K2,Image2] = Linear_Phase_Correction_Concatenate_new(R2,Rc2);
C1 = zeros(size(KK));
C2 = C1;
for j=-(size(C1,2))/2:size(C1,2)/2-1
    C1(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*0);
    C2(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*(mb-1));
end
slice_Grappa2