function [rawdata, ktraj_adc, adc_len, readouts, channels] = Read_RampSampled_EPI_Data(nF, path)
% path = 'E:\SMS-CAIPI\RampSampling_PC\revised_ramp_test_aug24\mb2';
% nF=6; % the number of the data set to load; you can grid any data with pulseq, but in this tutorial it is only necessary for #12
pattern='*.seq';
D=dir([path filesep pattern]);
[~,I]=sort(string({D(:).name}));
seq_file_path=[path filesep D(I(nF)).name];
% fprintf(['loading `' seq_file_path '´ ...\n']);
seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
[p,n,e] = fileparts(seq_file_path);
basic_file_path=fullfile(p,n);
data_file_path= [basic_file_path '.mat']; % try to load a matlab file with raw data...
try
    % fprintf(['loading `' data_file_path '´ ...\n']);
    data_unsorted = load(data_file_path);
    if isstruct(data_unsorted)
        fn=fieldnames(data_unsorted);
        assert(length(fn)==1); % we only expect a single variable
        data_unsorted=double(data_unsorted.(fn{1}));
    end
catch
    data_file_path= [basic_file_path '.dat']; % now try to load a raw data file...
    % fprintf(['falling back to `' data_file_path '´ ...\n']);
    twix_obj = mapVBVD(data_file_path);
    if iscell(twix_obj)
        data_unsorted = double(twix_obj{end}.image.unsorted());
        seqHash_twix=twix_obj{end}.hdr.Dicom.tSequenceVariant;
    else
        data_unsorted = double(twix_obj.image.unsorted());
        seqHash_twix=twix_obj.hdr.Dicom.tSequenceVariant;
    end    
    if length(seqHash_twix)==32
        % fprintf(['raw data contain pulseq-file signature ' seqHash_twix '\n']);
    end
    clear twix_obj
end

traj_recon_delay=0.2e-06;% [0.527 -1.367 0]; % adjust this parameter to potentially improve resolution & geometric accuracy. It can be calibrated by inverting the spiral revolution dimension and making two images match. for our Prisma and a particular trajectory we found 1.75e-6
grad_offsets=[0 0 0];

seq = mr.Sequence();              % Create a new sequence object
seq.read(seq_file_path,'detectRFuse');
seqName=seq.getDefinition('Name');
if ~isempty(seqName), fprintf('sequence name: %s\n',seqName); end
%[ktraj_adc, ktraj, t_excitation, t_refocusing, t_adc] = seq.calculateKspace('trajectory_delay', traj_recon_delay);
[ktraj_fov, ktraj_adc_fov, ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP('trajectory_delay', traj_recon_delay, 'gradient_offset', grad_offsets);
figure; plot(ktraj(1,:),ktraj(2,:),'b',...
             ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % a 2D plot
axis('equal'); 
title('2D k-space trajectory'); drawnow;
%% Define FOV and resolution and simple off-resonance frequency correction 

fov=256e-3; Nx=128; Ny=Nx; % define parameters explicitly
Ns=1; % this is a number of slices (or contrasts) which needs to be specified manually for now
%Na=1;

def_fov=seq.getDefinition('FOV'); % try to read from definitions 
if numel(def_fov)
    fov=max(def_fov);
end
%fov=fov*1.5;
deltak=1/fov;

% or estimate from the k-space trajectory
k_max=max(vecnorm(ktraj_adc));
Nx=round(k_max/deltak*2);
Ny=Nx; 

os=2; % oversampling factor (we oversample both in image and k-space)
offresonance=0; % global off-resonance in Hz

rawdata = permute(data_unsorted, [1,3,2]);
adc_len=size(rawdata,1);
readouts=size(rawdata,2);
channels=size(rawdata,3);
Na=numel(rawdata(:,:,1))/Ns/numel(t_adc); % averages/acquisitions
if Na>1
    nTRs=size(rawdata,2)/Na;
    rawdata = sum(permute(reshape(rawdata, [size(rawdata,1),size(rawdata,2)/Na,Na,size(rawdata,3)]),[1,2,4,3]),4);    
    readouts=readouts/Na;
end

if strcmp('ute_rs',seq.getDefinition('Name'))
    % average every 2nd spoke because of the half-rf excitation
    ktraj_adc=reshape(ktraj_adc,[3,adc_len,readouts]);
    ktraj_adc=ktraj_adc(:,:,1:2:end-1);
    t_adc=reshape(t_adc,[1,adc_len,readouts]);
    t_adc=t_adc(:,:,1:2:end-1);
    %rawdata=rawdata(:,1:2:end-1,:);
    rawdata=rawdata(:,1:2:end-1,:)+rawdata(:,2:2:end,:);
    readouts=readouts/2;
    ktraj_adc=reshape(ktraj_adc,[3,adc_len*readouts]);
    t_adc=reshape(t_adc,[1,adc_len*readouts]);
end

rawdata = reshape(rawdata, [size(rawdata,1)*size(rawdata,2)/Ns,Ns,size(rawdata,3)]);
ktraj_adc=ktraj_adc(:,1:end/Ns);
t_adc=t_adc(1:end/Ns);

for s=1:Ns
    for c=1:channels
        rawdata(:,s,c) = rawdata(:,s,c) .* exp(-1i*2*pi*t_adc'*offresonance);
    end
end
end