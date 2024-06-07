function image = f2i(epi_kxkyzc_lpcCor)

% image=sos(ifft2c(epi_kxkyzc_lpcCor(:,:,1,:)));

A=zeros(size(epi_kxkyzc_lpcCor));
for i=1:size(A,3)
    A(:,:,i)=fftshift(ifft2(ifftshift(epi_kxkyzc_lpcCor(:,:,i))));end
image=(sum(abs(A).^2,3)).^0.5;
figure();imshow(image,[])
end