function dc = epiphasecorrect(d, a)
% function dc = epiphasecorrect(d, a)
%
% Odd/even phase correction for EPI
%
% Inputs:
%   d      [nx, etl, ...]   Raw EPI data (Cartesian)
%   a      [2 1]            Odd/even phase mismatch linear fit parameters,
%                           obtained with getoephase.m using the same EPI train
%                           but with phase-encodes off.
%
% Outputs:
%   dc     size(d)          Input data after subtracting odd/even phase mismatch

nx = size(d,1);
etl = size(d,2);

dSize = size(d);

% Collapse higher dimensions for looping
d = reshape(d, nx, etl, []);

% We'll do correction in image space
x = fftshift(ifft(fftshift(d,1), [], 1),1);

% Construct 2D map of odd/even mismatch
X = [(-nx/2+0.5):(nx/2-0.5)]'/nx * ones(1, etl);  % see also getoephase.m
TH = a(1)*ones(size(X)) + a(2)*X;

if size(a,2) > 1
    % 'Secret' way to use this function:
    % use the measured odd/even phase mismatch instead of global linear fit
    TH = a;
else
    TH = TH(:,2:2:end);
end

% Subtract phase mismatch from even echoes
for ii = 1:size(d,3)
    x(:,2:2:end,ii) = x(:,2:2:end,ii).*exp(-1i*TH);
end

% Go back to data space and reshape
dc = fftshift(fft(fftshift(x,1), [], 1),1);
dc = reshape(dc, [nx  etl dSize(3:end)]);