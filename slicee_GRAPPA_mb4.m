function I = slicee_GRAPPA_mb4(KK,K1,K2,K3,K4)
mb=4;
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
end