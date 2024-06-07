function mask = new_Th(in,c,thereshold,degree)
mask = zeros(size(in,1),degree);
%amp = abs(in);
Nread = size(in,1)/c;
for j=1:c
    data = in((j-1)*Nread+1:j*Nread,:);
    % Max = data(:,1:3); Max = Max(:); Max = sort(Max,'descend'); Max = mean(Max(1:10));
    Mask1 = ones(Nread,1);
    % L_low = zeros(1,Nphase);
    % L_high = zeros(1,Nphase);
    for k=1:3
        Line = abs(data(:,k));
        % m = double(Line<=thereshold*Max);
        m = double(Line<=thereshold*max(Line));
        Mask1 = Mask1.*m;
    end
    m = smoothdata(Mask1,"gaussian",15);
    m = m>=0.9;
    mask((j-1)*Nread+1:j*Nread,:) = repmat(double(m),1,degree);
end
end
