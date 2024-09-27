%% set system limits
sys = mr.opts('MaxGrad', 28, 'GradUnit', 'mT/m', ...
    'MaxSlew', 160, 'SlewUnit', 'T/m/s', ... 
    'rfRingdownTime', 20e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 10e-6);

%% basic parameters
seq=mr.Sequence(sys);           % Create a new sequence object
fov = [256, 256, 140]*1e-3;
Nx = 128;
Ny = 128;
Nz = 64;
mb = 4;
alpha=70;                       % flip angle
sliceThickness = fov(3)/Nz;            % slice
slicesPerBand=Nz/mb;
bandStep = sliceThickness*slicesPerBand/2;
%TR=21e-3;                      % ignore TR, go as fast as possible
%TE=60e-3;                      % ignore TE, go as fast as possible
% more in-depth parameters
pe_enable=1;                    % a flag to quickly disable phase encoding (1/0) as needed for the delay calibration
rfDuration=8e-3;
roDuration=640e-6;              % not all values are possible, watch out for the checkTiming output
%% Create alpha-degree slice selection pulse and corresponding gradients 
[rf1, gz, gzReph] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,...
    'SliceThickness',sliceThickness,'apodization',0.42,'timeBwProduct',4,'system',sys);
tt=rf1.t;
rf1_center=mr.calcRfCenter(rf1);
signal_rf = zeros(size(rf1.signal));
for band=1:1
    % dF  = gz.amplitude*bandStep*(-(mb-1)/2+(band-1/2));
    dF  = gz.amplitude*bandStep*(band-1);
    % dPh = 2*pi/mb*band;
    dPh = 0;
    signal_rf = signal_rf+rf1.signal.*exp(1i*(2*pi*dF*(tt-rf1_center)+dPh));
end
rf = rf1;
rf.signal = signal_rf;
%% fat saturation
sat_ppm=-3.45;
sat_freq=sat_ppm*1e-6*sys.B0*sys.gamma;
rf_fs = mr.makeGaussPulse(110*pi/180,'system',sys,'Duration',8e-3,...
    'bandwidth',abs(sat_freq),'freqOffset',sat_freq,'use','saturation');
rf_fs.phaseOffset=-2*pi*rf_fs.freqOffset*mr.calcRfCenter(rf_fs); % compensate for the frequency-offset induced phase    
gz_fs = mr.makeTrapezoid('z',sys,'delay',mr.calcDuration(rf_fs),'Area',1/1e-4); % spoil up to 0.1mm

%% spoiler
spoiler_amp=3*8*42.58*10e2;
est_rise=500e-6; % ramp time 280 us
est_flat=2500e-6; %duration 600 us

gp_r=mr.makeTrapezoid('x','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys);
gp_p=mr.makeTrapezoid('y','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys);
gp_s=mr.makeTrapezoid('z','amplitude',spoiler_amp,'riseTime',est_rise,'flatTime',est_flat,'system',sys);

gn_r=mr.makeTrapezoid('x','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys);
gn_p=mr.makeTrapezoid('y','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys);
gn_s=mr.makeTrapezoid('z','amplitude',-spoiler_amp,'delay',mr.calcDuration(rf_fs), 'riseTime',est_rise,'flatTime',est_flat,'system',sys);

%%
[bw,f0,M_xy_sta,F1]=mr.calcRfBandwidth(rf);
[M_z,M_xy,F2]=mr.simRf(rf);
%
figure; plot(F1,abs(M_xy_sta),F2,abs(M_xy),F2,M_z);
axis([f0-2*bw, f0+2*bw, -0.1, 1.2]);
legend({'M_x_ySTA','M_x_ySIM','M_zSIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('STA vs. simulation, flip angle 70°');

figure; plot(F2,atan2(abs(M_xy),M_z)/pi*180);
axis([f0-2*bw, f0+2*bw, -5, 100]);
xlabel('frequency offset / Hz');
ylabel('flip ange [°]');
legend({'SINC'});
grid on;
title('Achieved flip angle for the nominal 70° flip');

figure; plot(F2,real(M_xy),F2,imag(M_xy));
axis([f0-2*bw, f0+2*bw, -1.2, 1.2]);
legend({'M_xSIM','M_ySIM'});
xlabel('frequency offset / Hz');
ylabel('magnetisation');
title('Real and imag. parts of transverse magnetisation, 70° flip');
%%
% Define other gradients and ADC events
deltak=1./fov; % Pulseq default units for k-space are inverse meters
kWidth = Nx*deltak(1);
Ablip = 1/bandStep/mb;
A = max((mb-1)*Ablip,deltak(2));
% blip_dur = ceil(2*sqrt(deltak/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2;
% start with the blip
% blip_dur = 1.5*2*ceil(2*sqrt(A/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime;
blip_dur = 2*ceil(2*sqrt(A/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime;
gyBlip = mr.makeTrapezoid('y',sys,'Area',-deltak(2),'Duration',blip_dur); % we use negative blips to save one k-space line on our way towards the k-space center
gzBlip = mr.makeTrapezoid('z',sys,'Area',(mb-1)*Ablip,'Duration',blip_dur);
%gzBlip = mr.makeTrapezoid('z',sys,'Area',-(mb-1)*deltak(3),'Duration',blip_dur);

% blip_dur = ceil(2*sqrt((mb-1)*deltak(3)*bandStep/sliceThickness/sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2;
% blip_dur = ceil(2*sqrt((mb-1)*deltak(3)*sys.maxSlew)/sys.gradRasterTime/2)*sys.gradRasterTime*2;

extra_area=blip_dur/2*blip_dur/2*sys.maxSlew;
gx = mr.makeTrapezoid('x',sys,'Area',kWidth+extra_area,'duration',roDuration+blip_dur);
actual_area=gx.area-gx.amplitude/gx.riseTime*blip_dur/2*blip_dur/2/2-gx.amplitude/gx.fallTime*blip_dur/2*blip_dur/2/2;
gx=mr.scaleGrad(gx,kWidth/actual_area);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'system',sys); % if no 'Duration' is provided shortest possible duration will be used
gyPre = mr.makeTrapezoid('y','Area',(Ny/2-1)*deltak(2),'system',sys);
 
adcDwellNyquist=deltak(1)/gx.amplitude; % dwell time on the top of the plato
adcDwell=floor(adcDwellNyquist/sys.adcRasterTime)*sys.adcRasterTime;
adcSamples=floor(roDuration/adcDwell/4)*4; % on Siemens the number of ADC samples need to be divisible by 4
adc = mr.makeAdc(adcSamples,'Dwell',adcDwell);

time_to_center=adc.dwell*((adcSamples-1)/2+0.5); % Pulseq (and Siemens) define the samples to happen in the center of the dwell period
adc.delay=round((gx.riseTime+gx.flatTime/2-time_to_center)/sys.rfRasterTime)*sys.rfRasterTime; 

% split the blip into two halves and produce a combined synthetic gradient
gyBlip_parts = mr.splitGradientAt(gyBlip, blip_dur/2, sys);
[gyBlip_up,gyBlip_down,~]=mr.align('right',gyBlip_parts(1),'left',gyBlip_parts(2),gx);
% now for inner echos create a special gy gradient, that will ramp down to 0, stay at 0 for a while and ramp up again
gyBlip_down_up=mr.addGradients({gyBlip_down, gyBlip_up}, sys);

% split the blip into two halves and produce a combined synthetic gradient
gzBlip_parts = mr.splitGradientAt(gzBlip, blip_dur/2, sys);
[gzBlip_up,gzBlip_down,~]=mr.align('right',gzBlip_parts(1),'left',gzBlip_parts(2),gx);
% now for inner echos create a special gy gradient, that will ramp down to 0, stay at 0 for a while and ramp up again
gzBlip_down_up=mr.addGradients({mr.scaleGrad(gzBlip_down,1/(mb-1)), mr.scaleGrad(gzBlip_up,1/(mb-1))}, sys);
gzBlip_down_down=mr.addGradients({mr.scaleGrad(gzBlip_down,1/(mb-1)), mr.scaleGrad(gzBlip_up,-1)}, sys);
gzBlip_up_up=mr.addGradients({mr.scaleGrad(gzBlip_down,-1), mr.scaleGrad(gzBlip_up,1/(mb-1))}, sys);

% pe_enable support
gyBlip_up=mr.scaleGrad(gyBlip_up,pe_enable);
gyBlip_down=mr.scaleGrad(gyBlip_down,pe_enable);
gyBlip_down_up=mr.scaleGrad(gyBlip_down_up,pe_enable);
gyPre=mr.scaleGrad(gyPre,pe_enable);

% gradient spoiling
gzSpoil=mr.makeTrapezoid('z','Area',4/sliceThickness,'system',sys); % 4 cycles over the slice thickness

% skip timing (TE/TR calculation), we'll accept the shortest TE/TR
%gzPre1 = mr.makeTrapezoid('z','Area',-3*deltak(3)*bandStep/sliceThickness,'system',sys);
% define sequence blocks

seq.addBlock(gp_r,gp_p,gp_s);
seq.addBlock(rf_fs, gn_r,gn_p,gn_s);
seq.addBlock(rf,gz);
% gzReph = mr.makeTrapezoid('z','Area',gzReph.area-gzBlip.area/2,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',gzReph.area,'system',sys);
seq.addBlock(mr.align('left',gyPre,gzReph,'right',gxPre));
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();
%A = ktraj(3,end);
% gz_compensate = mr.makeTrapezoid('z',sys,'Area',-A,'Duration',blip_dur);
% seq.addBlock(gz_compensate);
% seq.addBlock(mr.scaleGrad(gzBlip,-1/2));
%seq.addBlock(mr.scaleGrad(gzBlip,-1/2),mr.makeDelay(blip_dur));
%seq.addBlock(mr.makeDelay(blip_dur));
%seq.addBlock(mr.makeDelay(mr.calcDuration(gzPre1)));
for i=1:Ny % loop over phase encodes
    if i==1
        seq.addBlock(gx,gyBlip_up,adc); % Read the first line of k-space with a single half-blip at the end
    elseif i==Ny
        seq.addBlock(gx,gyBlip_down,adc); % Read the last line of k-space with a single half-blip at the beginning
    elseif (mod(i,mb)==1 && i~=1)
        seq.addBlock(gx,gyBlip_down_up,adc);
    elseif (mod(i,mb)==0 && i~=Ny)
        seq.addBlock(gx,gyBlip_down_up,adc);
    elseif (mod(i,mb)~=0 && mod(i,mb)~=1)
        seq.addBlock(gx,gyBlip_down_up,adc);
    end 
    gx = mr.scaleGrad(gx,-1);   % Reverse polarity of read gradient),'Duration',mr.calcDuration(gxPre),'system',sys);
end
seq.addBlock(gzSpoil);

%% check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;

if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% prepare sequence export
seq.setDefinition('FOV', [fov fov sliceThickness]);
seq.setDefinition('Name', 'S1');
seq.write('S1.seq')       % Write to pulseq file
%seq.install('siemens');

%% plot sequence and k-space diagrams

%seq.plot();
seq.plot('timeDisp','us','showBlocks',1); %detailed view

% k-space trajectory calculation
[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = seq.calculateKspacePP();

% plot k-spaces
figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D k-space plot
axis('equal'); % enforce aspect ratio for the correct trajectory display
hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
title('full k-space trajectory (k_x x k_y)');
figure();plot3(ktraj_adc(1,:),ktraj_adc(2,:),ktraj_adc(3,:),'-k')

%% PNS calc

[pns_ok, pns_n, pns_c, tpns]=seq.calcPNS('idea/asc/MP_GPA_K2309_2250V_951A_AS82.asc'); % prisma
if (pns_ok)
    fprintf('PNS check passed successfully\n');
else
    fprintf('PNS check failed! The sequence will probably be stopped by the Gradient Watchdog\n');
end

%% very optional slow step, but useful for testing during development e.g. for the real TE, TR or for staying within slewrate limits  
rep = seq.testReport;
fprintf([rep{:}]);
