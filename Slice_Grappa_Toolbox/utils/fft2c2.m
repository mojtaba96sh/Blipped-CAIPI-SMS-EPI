function res = fft2c2(x)
fctr = size(x,1)*size(x,2);

res = zeros(size(x));

for m = 1:size(x,4)
    for n=1:size(x,3)
        res(:,:,n,m) = 1/sqrt(fctr)*fftshift(fft2(ifftshift(x(:,:,n,m))));
    end
end


