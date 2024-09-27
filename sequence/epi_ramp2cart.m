function [cartData] = epi_ramp2cart(rampData,ktraj_adc,Ny)

rampData = squeeze(rampData);
Ny_prim = Ny;
Nx_prim = size(ktraj_adc,2)/Ny;
C = size(rampData,2);
Nx = Ny;
cartData = zeros(Nx,C,Ny);
% ktraj_adc_new = zeros();
for i=1:Ny
    for c=1:C
        K_ramp = rampData((i-1)*Nx_prim+1:i*Nx_prim,c);
        kx_ramp = ktraj_adc(1,(i-1)*Nx_prim+1:i*Nx_prim);
        % ky_ramp = ktraj_adc(2,(i-1)*Nx_prim+1:i*Nx_prim);

        kx_cart = linspace(min(kx_ramp),max(kx_ramp),Nx);
        % ky_cart = ky_ramp(Ny/2)*ones(1,Ny);

        % K_cart = interp1(kx_ramp, K_ramp, kx_cart, 'spline');
        K_cart = interp1(kx_ramp, K_ramp, kx_cart, 'linear');
        % K_cart = griddata(kx_ramp, ky_ramp, K_ramp, kx_cart, ky_cart);
        cartData(:,c,i) = K_cart;
    end
end
end
