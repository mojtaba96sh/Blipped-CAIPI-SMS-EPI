function I = slicee_GRAPPA_recon(KK,K1,K2)
mb=2;
C1 = zeros(size(KK));
C2 = C1;
for j=-(size(C1,2))/2:size(C1,2)/2-1
    C1(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*0);
    C2(:,j+size(C1,2)/2+1,:) = exp(-1i*2*pi*j/mb*(mb-1));
end
slice_Grappa2
end