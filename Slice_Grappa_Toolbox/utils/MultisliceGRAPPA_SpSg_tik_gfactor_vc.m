function [K, w, G_factor] = MultisliceGRAPPA_SpSg_tik_gfactor_vc(K_Collapsed, varargin)

%K = MultisliceGRAPPA....(K_Collapsed,w,KernelSize)
%K = MultisliceGRAPPA....(K_Collapsed,K_Indiv,KernelSize,DataForFit,prot)

% October 24 2009
% Kawin Setsompop

% add leakage block (sp-SG) kernel calc to reduce leakage artifact
% Kawin Setsompop
% June 2012

% Inputs:
% K_Collapsed  : the slice collapsed k space data 
% K_Indiv      : the 'acs' data. the data from individually acquired slice (for one simultanous multislice group)
% KernelSize   : [ NumberOfElementsAlongKx NumberOfElementsAlongKy distanceBetweenElementsAlongKx distanceBetweenElementsAlongdKy]
% DataForFit   : Part of the acs data that will be use in grappa kernel generation 
%                This can be specified in 2 ways:
%                DataForFit = 'full'; use all the data from the acs dataset to estimate the grappa kernel
%                DataForFit = [SizeAlongKx SizeAlongKy]; use the center part of the acs k-space of size [SizeAlongKx SizeAlongKy] to estimate the kernel
% prot         : an array containing protocol parameters. In this function only  'prot.lPhaseEncodingLines' is used. This represents the Number of phase encoding line in the acquisition (including the missing line due to P.F.)

% Outputs:
% K            : the slice seperated dataset
% w            : kernels 
% 
% K space data matrix dimensions are based on standard siemens convention
%   [#01]  NColMeas: 
%   [#02]  NLinMeas: 
%   [#03]  NChaMeas: 
%   [#04]  NSetMeas: 
%   [#05]  NEcoMeas: 
%   [#06]  NPhsMeas: 
%   [#07]  NRepMeas: 
%   [#08]  NSegMeas: 
%   [#09]  NParMeas: 
%   [#10]  NSlcMeas: 
%   [#11]  NIdaMeas: 
%   [#12]  NIdbMeas: 
%   [#13]  NIdcMeas: 
%   [#14]  NIddMeas: 
%   [#15]  NIdeMeas: 
%   [#16]  NAveMeas:   



Nchannels = size(K_Collapsed,3);

if size(varargin,2) <= 3
    % apply kernel in recon
    
    w = varargin{1};
    KernelSize = varargin{2};
    KernelSupplied = 1;
    NslicesEX = size(w,3);

    if size(varargin,2) == 3
        virtual = varargin{3};
    else
        virtual = 0;
    end
    
else
    % estimate kernel and g-factor
    
    KernelSupplied = 0;
    K_Indiv = varargin{1};
    KernelSize = varargin{2};
    DataForFit = varargin{3};    
    prot = varargin{4};
    
    lambda_tik = varargin{5}
    
    Weights_p = varargin{6};
    
    if size(varargin,2) == 7
        virtual = varargin{7};
    else
        virtual = 0;
    end

    
    NslicesEX = max(size(K_Indiv));
end

if length(KernelSize) == 2 % dim 3 and 4 are for SkipSteps
    KernelSize(3) = 1;
    KernelSize(4) = 1;
end

if KernelSupplied == 0
    
    disp(['Estimating kernel'])
    
    s = size(K_Indiv{1});
    K_CollapsedFake = zeros(s);
    for count = 1:NslicesEX
        K_CollapsedFake = K_CollapsedFake + K_Indiv{count};
    end
    
    % create Indexing
    if isfloat(DataForFit) % use only part of the supplied ACS data
        Nread = size(K_Collapsed,1);
        Nlin_full = prot.lPhaseEncodingLines;
        Nlin_part = size(K_Collapsed,2);
        
        DataEdgeBegin = [1 1] + floor(([Nread Nlin_full] - DataForFit)/2);
        DataEdgeEnd   = [Nread Nlin_full] - ceil(([Nread Nlin_full] - DataForFit)/2);
        
        if Nlin_full ~= Nlin_part % i.e. Partial Fourier so need to make sure use data in the center area (around dc) for kernel fitting
            pad_lines = Nlin_full - Nlin_part;
            DataEdgeBegin(2) = DataEdgeBegin(2) - pad_lines; % assuming missing part is 'pre'
            DataEdgeEnd(2) = DataEdgeEnd(2) - pad_lines;
            if DataEdgeBegin(2) < 1
                display('fitting data area exceed acquired data because of P.F')
                keyboard
            end
        end
    else
        % use full ACS data
        DataForFit = [size(K_Collapsed,1) size(K_Collapsed,2)];
        DataEdgeBegin = [1 1] ;
        DataEdgeEnd   = DataForFit;
    end
    
    
    KxIndex = -floor( KernelSize(1)/2):floor((KernelSize(1)-1)/2); % i.e. kernel will be bulking in the "up-left" direction if kernel_size is even
    KyIndex = -floor( KernelSize(2)/2):floor((KernelSize(2)-1)/2);
    KxIndex = KxIndex*KernelSize(3);
    KyIndex = KyIndex*KernelSize(4);
        
    BeginIndex = DataEdgeBegin  + floor(KernelSize(1:2)/2).*[KernelSize(3) KernelSize(4)]; %current mid point of the kernel
    EndIndex = DataEdgeEnd - floor((KernelSize(1:2)-1)/2).*[KernelSize(3) KernelSize(4)];
    
    % Fitting Kernel
    
    %Calculate DataMatrix
    tic
    
    MatrixSizeSliceSeperation = [(DataForFit(1) - (KernelSize(1)-1)*KernelSize(3))*(DataForFit(2) - (KernelSize(2)-1)*KernelSize(4)) ,     ( KernelSize(1) * KernelSize(2) )*Nchannels ]

    % tell you how well conditioned the problem is i.e. get good Kernel weights if m>>n where [m,n] = size(DataMatrix)
    
    
    DataMatrix = zeros( MatrixSizeSliceSeperation .* [NslicesEX,1] );
    
    
    for SliceCount = 1:NslicesEX
        DataMatrixCurrent = zeros(MatrixSizeSliceSeperation);
        
        CurrentIndex = BeginIndex;
        count = 1;
        
        for countY = 1:(DataForFit(2) - (KernelSize(2)-1)*KernelSize(4))
            for countX = 1:(DataForFit(1) - (KernelSize(1)-1)*KernelSize(3))
                DataMatrixCurrent(count,:) = reshape( K_Indiv{SliceCount}(KxIndex+CurrentIndex(1),KyIndex+CurrentIndex(2),:) ,1,[],1);
                count = count+1;
                CurrentIndex = CurrentIndex + [1 0];
            end
            CurrentIndex(1) = BeginIndex(1);
            CurrentIndex(2) = CurrentIndex(2) + 1;
        end
         
        
        DataMatrix( 1 + (SliceCount-1) * MatrixSizeSliceSeperation(1) : SliceCount *MatrixSizeSliceSeperation(1), : ) = DataMatrixCurrent;
        
    end
       
    
    disp(['DataMatrix Form: ' num2str(toc) ' s'])
    
    % Calculate Weights
    tic
    
    if ~lambda_tik
        InvDataMatrix = pinv(DataMatrix);
    else    
        [u,s,v] = svd(DataMatrix, 'econ');

        s_inv = diag(s);
        
        disp(['Min sigma: ', num2str(min(s_inv)), '   Max sigma: ', num2str(max(s_inv))])
        
        s_inv = conj(s_inv) ./ (abs(s_inv).^2 + lambda_tik);

        InvDataMatrix = v * diag(s_inv) * u'; 
    end
    
    disp(['Inversion: ' num2str(toc) ' s'])
    

    % kernel to recon only actual coils
    w = zeros([Nchannels * prod(KernelSize), Nchannels / (1+virtual), NslicesEX]);
    
    flag = 1;
    for ChCount = 1 : (Nchannels / (1+virtual))
        for SliceCount = 1:NslicesEX
            
            ACS_actual = K_Indiv{SliceCount}(BeginIndex(1):EndIndex(1),BeginIndex(2):EndIndex(2),ChCount);
            
            if flag == 1
                sV = size(ACS_actual(:));
                flag = 0;
            end
            
            ACS = cat(1,zeros(sV(1)*(SliceCount-1),1),ACS_actual(:), zeros(sV(1)*(NslicesEX -SliceCount),1));
            
            w(:,ChCount,SliceCount) = InvDataMatrix*ACS(:);
        end
    end
    
    
end



K = 0;



if KernelSupplied~=0
    disp('Applying kernel')


    % Seperate Slices using the calculated Kernel
    s = size(K_Collapsed); 
    if length(s) < 10
        s(end+1:9) = 1;
    end
    s(10) = NslicesEX;   
    K = zeros(s);

    Wcurrent = zeros((KernelSize(1)-1)*KernelSize(3)+1, (KernelSize(2)-1)*KernelSize(4)+1, Nchannels);

    % recon only actual channels
    
    disp(['Channels to recon: ', num2str(Nchannels / (1+virtual))])
    
    %tic
    for SliceCount = 1:NslicesEX
        disp(['Recon slice: ', num2str(SliceCount), ' / ', num2str(NslicesEX)])
        
        for ChCount = 1 : (Nchannels / (1+virtual))
            
            Wcurrent(1:KernelSize(3):end,1:KernelSize(4):end,:) = reshape(w(:,ChCount,SliceCount), KernelSize(1),KernelSize(2),Nchannels);    
            
            for SetCount = 1:size(K_Collapsed,4)
                for AvgCount = 1:size(K_Collapsed,7)
                    
                    for ChCount2 = 1:Nchannels 
                    
                        WcurrentCh2 = Wcurrent(end:-1:1,end:-1:1,ChCount2);
                        
                        K(:,:,ChCount,SetCount,:,:,AvgCount,:,:,SliceCount) =  K(:,:,ChCount,SetCount,:,:,AvgCount,:,:,SliceCount) + ...
                                     conv2(K_Collapsed(:,:,ChCount2,SetCount,:,:,AvgCount,:,:,:), WcurrentCh2,'same') ;
                                 
                    end
                end
            end
        end
    end

end







if nargout > 2 
    % image space grappa weights

    kernel_hsize = (KernelSize-1)/2;

    N = [Nread, Nlin_full];
    
    G_factor = zeros([N, NslicesEX]);

    % number of actual coils:
    num_actual = size(Weights_p, 3);

    
    for slc = 1:NslicesEX
    
        disp(['G-factor slice: ', num2str(slc), ' / ', num2str(NslicesEX)])
        
        w_slice = w(:,:,slc);

        image_weights = zeros([N, Nchannels, num_actual], 'single');
 

        for c = 1:num_actual

            image_weights(1+end/2, 1+end/2, c, c) = 1;

            weights = w_slice(:, c);

            W = reshape(weights, [KernelSize, Nchannels]);

            for coil = 1:Nchannels
                image_weights( 1+end/2 - kernel_hsize(1) : 1+end/2 + kernel_hsize(1), 1+end/2 - kernel_hsize(2) :  1+end/2 + kernel_hsize(2), coil, c) = W(:,:,coil);
            end

        end

        image_weights = flipdim( flipdim( ifft2c2( image_weights ), 1 ), 2 ) * sqrt(prod(N));
        

        weights_p = Weights_p(:,:,:,slc);


        g_cmb = zeros([N, Nchannels]);

        for c = 1:Nchannels
            weights_coil = squeeze( image_weights(:,:,c,1:num_actual) );

            g_cmb(:,:,c) = sum(weights_p .* weights_coil, 3);
        end

        g_comb = sqrt( sum(g_cmb .* conj( g_cmb ), 3) );

        p_comb = sqrt( sum(weights_p .* conj( weights_p ), 3) );

        G_factor(:,:,slc) = g_comb ./ p_comb;
        
    end
    
end



end
