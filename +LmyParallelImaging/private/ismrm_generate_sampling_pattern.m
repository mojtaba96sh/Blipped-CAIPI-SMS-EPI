function pat = ismrm_generate_sampling_pattern(data_shape, acc, ref, sshift)
%
%  [data, sampling_pattern] = ismrm_sample_data(img_obj, csm, acc, ref)
%
%  Samples the k-space of object provided in 'img_obj' after first applying
%  coil sensitivity maps in 'csm' and Fourier transforming to k-space.
%
%  INPUT:
%    - img_obj [x,y]    : Object in image space
%    - csm     [x,y,c]  : Coil sensitivity maps
%    - acc     scalar   : Acceleration factor
%    - ref     scalar   : Reference lines (in center of k-space)
%    - sshift  scalar   : Sampling shift, i.e for undersampling, do we
%                         start with line 1 or line 1+sshift?
%
%  OUPUT:
%    - data    [kx,ky,c]: Sampled data in k-space (zeros where not sampled)
%    - pat     [kx,ky,c]: Sampling pattern (0 = not sampled,
%                                           1 = imaging data,
%                                           2 = reference data,
%                                           3 = reference and imaging data)
%
%
%   Code made available for the ISMRM 2013 Sunrise Educational Course
% 
%   Michael S. Hansen (michael.hansen@nih.gov)
%

if nargin < 1,
    error('At least two arguments are needed');
end

if nargin < 2,
    acc = 1;
end

if nargin < 3,
    ref = 0;
end

if nargin < 4,
    sshift = 0;
end

sshift = mod(sshift,acc);


%Generate parallel imaging undersampling pattern
pat_img = zeros(data_shape);
pat_img(:,(1+sshift):acc:end) = 1;

%Generate reference lines pattern
pat_ref = zeros(data_shape);
if (ref > 0),
    pat_ref(:,(1:ref)+bitshift(data_shape(2),-1)-bitshift(ref,-1)) = 2;
end

pat = pat_img + pat_ref;

%Fully sampled data in kspace:
%data = ismrm_transform_image_to_kspace(repmat(img_obj, [1 1 size(csm,3)]) .* csm, [1 2]);

%Apply sampling pattern:
%data = data .* repmat((pat > 0), [1 1 size(csm,3)]); 

return 