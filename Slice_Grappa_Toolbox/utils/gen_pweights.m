function [ weights_p ] = gen_pweights( kspace_acs, acs_size )
%GEN_PWEIGHTS Summary of this function goes here
%   Detailed explanation goes here


% Tukey window on acs image to create rsos coil combination weights

tw1 = tukeywin(acs_size(1));
tw2 = tukeywin(acs_size(2));

[N(1), N(2), num_chan] = size(kspace_acs);


tw = repmat( padarray(tw1 * tw2', (N-acs_size)/2), [1,1,num_chan]);


img_acs = ifft2c( kspace_acs .* tw );


weights_p = conj(img_acs) ./ (eps + repmat(rsos(img_acs, 3), [1,1,num_chan]));



end

