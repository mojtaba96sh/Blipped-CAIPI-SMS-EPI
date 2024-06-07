%--------------------------------------------------------------------------
%% load data
%--------------------------------------------------------------------------

addpath utils/

load IMG

 
[N(1), N(2), num_chan, MB_factor, ] = size(IMG);


mosaic( squeeze( rsos(IMG, 3) ), MB_factor, 1, 1, '', genCaxis(rsos(IMG, 3)) *.8)

 

 
%%------------------------------------------------------------------------- 
%% collapse slices: 
%--------------------------------------------------------------------------

FOV_shift = 3;
    
shift_amount = round( (FOV_shift~=0) * [-(MB_factor/2 - 1) : MB_factor/2] * N(2) / (eps + FOV_shift) );
    
disp(num2str(shift_amount))



IMG_shift = zeros(size(IMG));

for s = 1:MB_factor
    IMG_shift(:,:,:,s) = circshift( IMG(:,:,:,s,:), [0,shift_amount(s),0,0] );
end

mosaic( squeeze( rsos(IMG_shift, 3) ), MB_factor, 1, 2, '', genCaxis(rsos(IMG, 3)) *.8)

%%------------------------------------------------------------------------- 
%% split slice grappa: g-factor
%--------------------------------------------------------------------------

lambda_tik = eps;               % tikhonov for kernel estimation

KernelSize = [7,7];         

DataForFit = [N(1),24];         % ACS size

prot.lPhaseEncodingLines = N(2);


% collapsed k-space
K_Collapsed = fft2c2( sum(IMG_shift, 4) );


% individual k-space
K_Indiv = [];

for s = 1:MB_factor
    K_Indiv{s} = fft2c2( IMG_shift(:,:,:,s) );
end



% coil combination weights for g-factor estimation 
Weights_p = zeros([N, num_chan, MB_factor]);

msk_acs = zeros([N, num_chan]);
msk_acs(1+end/2 - DataForFit(1)/2 : end/2 + DataForFit(1)/2, 1+end/2 - DataForFit(2)/2 : end/2 + DataForFit(2)/2, :) = 1;

for s = 1:MB_factor
    kspace_acs = K_Indiv{s} .* msk_acs;

    Weights_p(:,:,:,s) = gen_pweights( kspace_acs, DataForFit );
end


% g-factor and kernel estimation
[~, w, G_factor] = MultisliceGRAPPA_SpSg_tik_gfactor_vc(K_Collapsed, K_Indiv, KernelSize, DataForFit, prot, lambda_tik, Weights_p);



% apply kernel in k-space
K = MultisliceGRAPPA_SpSg_tik_gfactor_vc(K_Collapsed, w, KernelSize);    

IMG_recon = ifft2c2( squeeze( K ) );


% shift slices back
for s = 1:MB_factor
    IMG_recon(:,:,:,s) = circshift( IMG_recon(:,:,:,s), [0,-shift_amount(s),0] );

    G_factor(:,:,s) = circshift( G_factor(:,:,s), [0,-shift_amount(s)] );
end
    


   
rmse_sg = 100 * norm2( rsos(IMG_recon,3 ) - rsos(IMG,3) ) / norm2( rsos(IMG,3) );


mosaic( squeeze(rsos(IMG, 3)), 2, MB_factor/2, 50, ['R=1 reference data'], [0,0.5] )
mosaic( squeeze(rsos(IMG_recon, 3)), 2, MB_factor/2, 51, ['slice grappa rmse: ', num2str(rmse_sg), '%'], [0,0.5] )

mosaic( 1 ./ (eps + G_factor), 2, MB_factor/2, 101, ['max g-factor:', num2str(max(G_factor(:)))], [0,1] ), colormap jet

  
     
  