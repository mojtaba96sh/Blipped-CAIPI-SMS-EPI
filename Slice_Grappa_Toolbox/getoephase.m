function [a, th] = getoephase(x, verbose, threshold)
% function a = getoephase(x, [verbose=false], [threshold=0.1])
%
% Estimate odd/even EPI echo phase difference (linear and constant term)
% from one EPI echo train (without phase-encoding blips).
%
% Inputs
%   x    [nx etl nCoils]   EPI echo train of length 'etl', without phase-encoding.
%                          Assumes x has been inverse FFT'd along 1st (readout, or 'x') dimension,
%                          i.e., that x contains multiple 1D spatial profiles.
%                          Also assumes that etl is even.
%   verbose  boolean       Default: false
%
% Output
%   a    [2 1]     Linear fit to odd/even phase offsets
%                  a(1): constant phase offset (radians)
%                  a(2): linear term (radians/fov)
%                  The corresponding k-space shift (in samples) is a(2)/(2*pi)

if nargin < 2
    verbose = false;
end
if nargin < 3
    threshold = 0.1;   % magnitude threshold (max amplitude normalized to 1)
end

[nx etl nCoils] = size(x);

if mod(etl,2)
    error('etl must be even');
end

% Estimate off-resonance phase from even echoes (center line).
% Magnitude-squared coil weighting as in phase-contrast MRI.
th = zeros(1, etl/2);
for ic = 1:nCoils
    xe = x(end/2,2:2:etl,ic);  % even echoes
    tmp = unwrap(angle(xe));
    tmp = tmp - tmp(floor(end/2)); % need common (arbitrary) reference since we're combining coils
    th = th + abs(xe).^2 .* exp(1i*tmp);
end

if etl/2 < 10 
    maSpan = 1;
else
    maSpan = 5;   % moving-average span for smoothing
end

th = smooth(unwrap(angle(th)), maSpan);

% Subtract off-resonance phase from all echoes
if true
    % Fit straight line
    X = [2:2:etl]';
    B = [ones(length(th),1) X];
    a = B\th(:);  % a(2) = off-resonance phase accrual per echo (radians)
    XY = ones(nx,1) * [1:etl] ;  % [nx etl]
    DPH = a(1) + a(2)*repmat(XY, [1 1 nCoils]);  % linear fit to phase accrual (radians)
else
    % use measured off-resonance phase directly
    th = interp1(2:2:etl, th, 1:etl, 'linear', 'extrap');
    DPH = repmat(ones(nx,1) * th, [1 1 nCoils]);
end

%DPH = DPH - th(1);  % reference to first echo (arbitrarily)

xc = x.*exp(-1i*DPH);  % 'c' for 'corrected'

if verbose
    coil = ceil(nCoils/2);
    figure; subplot(121); im(angle(x(:,:,coil)));
    subplot(122); im(angle(xc(:,:,coil)));
    figure;
    plot(angle(x(end/2,:,coil)),'r'); hold on
    plot(angle(xc(end/2,:,coil)), 'g'); hold on
end

% spatial mask
% exclude object outside center FOV/2 (in x)
rssim = sqrt(sum(abs(xc).^2, 3));
mask = rssim > threshold*max(rssim(:));
mask([1:round(nx/4) round(3/4*nx):end],:) = false;

% Get odd/even phase mismatch for all neighboring echo pairs
th = zeros(nx, etl/2);
for ic = 1:nCoils
    xo = xc(:, 1:2:etl, ic);  % odd echoes
    xe = xc(:, 2:2:etl, ic);  % even echoes
    th = th + abs(xe).^2 .* exp(1i*angle(xe./xo)); % NB! We assume no phase-wrap
end
th = angle(th);  
if verbose
    figure; im(th.*mask(:,2:2:end), pi*[-1 1]); colormap default; colorbar;
    xlabel('position (pixel) along x'); ylabel('odd/even phase mismatch');
    title(sprintf('odd/even phase mismatch for all echo pairs (etl = %d)', etl));
end

% Do linear fit to the later (more stable) echo pairs
mask = mask(:, (end/2+1):etl);  % = size(th)
mask(:, 1:floor(end/2)) = false;
X = [(-nx/2+0.5):(nx/2-0.5)]'/nx * ones(1, size(th,2));
H = [ones(sum(mask(:)),1) X(mask)];  % spatial basis matrix
a = H\th(mask); 
if verbose
    % display range
    R = [min(th(mask)) max(th(mask))];

    figure; 

    subplot(221); im(th.*mask, R); colormap default; 
    title('measured odd/even phase difference'); colorbar;
    thhat = embed(H*a, mask);

    subplot(222); im(thhat.*mask, R); colormap default; 
    title('linear fit'); colorbar;

    subplot(223); im((thhat-th).*mask, 0.1*[-1 1]); colormap default; 
    title('difference'); colorbar;

    subplot(224);
    X = [(-nx/2+0.5):(nx/2-0.5)]'/nx * ones(1, etl);
    TH = a(1)*ones(size(X)) + a(2)*X;
    im(TH(:,(end/2+1):end).*mask, R); colormap default;
    title('a(1) + a(2)*x'); colorbar;
end

return