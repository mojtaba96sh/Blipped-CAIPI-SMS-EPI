function mask = Threshold(in,c,thereshold)
mask = zeros(size(in,1),2);
amp = abs(in);
Nread = size(in,1)/c;
for j=1:c
    data = amp((j-1)*Nread+1:j*Nread,:);
    Th = max(data(:))*thereshold;
    m = data>=Th;
    mask1 = double(sum(m,2)>0);
    mask((j-1)*Nread+1:j*Nread,:) = repmat(mask1,1,2);
end
end