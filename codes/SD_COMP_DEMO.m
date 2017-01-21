function [] = SD_COMP_DEMO(opts)
% This is a demo showing compressive sensing procedure can be inherently 
% incorporated into delta-sigma ADC. For the detailed explanation, please 
% consult with github readme file.
% In this script, one can choose different options which will be explained
% below:

if ~isempty(opts) && ~isstruct(opts)
    error('"opts" must be a structure');
end

function out = setOpts( field, default )
    if ~isfield(opts,field)
        opts.(field) = default;
    end
    out = opts.(field);
end
% Part I: Create and Save PHI matrix:
% With a given input, the random PHI matrix uniquely determined the
% compressive sensing results. To reproduce the results, one can save and
% load the PHI matrix. Also the quality of randomness of PHI should be
% checked to ensure reconstruction algorithm working.

clc;

PATHS = setOpts('PATHS',15);
DATA_LENGTH = setOpts('DATA_LENGTH',2^8);
WARMUP_CYCLES = setOpts('WARMUP_CYCLES',10);
SCALE = setOpts('SCALE',1);

numPointsFFT = DATA_LENGTH;

BITS = setOpts('BITS',1);

PHI = 2*(randn(PATHS,DATA_LENGTH)>0)-1;


% FS = setOpts('FS',10e6);
% TR = setOpts('TR',100e-12);
% PHI_wiggle = ((repmat(PHI,1,WARMUP_CYCLES)+1)/2)';
% PHI_stich = zeros(2*WARMUP_CYCLES*DATA_LENGTH,PATHS);
% TIME = (0:WARMUP_CYCLES*DATA_LENGTH)'/FS;

FILENAME = setOpts('FILENAME','TEST');

save_flag = setOpts('SAVE',1);
load_flag = setOpts('LOAD',0);

if save_flag==1
    save(FILENAME,'PHI');
end

if load_flag==1
    PHI = load(FILENAME,'PHI');
end
% plot(TIME_stich,PHI_stich)

binFFT = setOpts('binFFT',[3 51]);
amp = setOpts('amp',[0.8, 0.3]);
if length(amp)~=length(binFFT)
    error('amp and binFFT should have the same length. You beed to assign each frequency component an amplitude');
end

SIGNAL = zeros(1,numPointsFFT);
for ii = 1:length(binFFT)
    SIGNAL = SIGNAL+...
    (2^BITS/2*amp(ii)-0.005)*sin(2*pi*(0:numPointsFFT-1)/numPointsFFT*binFFT(ii));
end
SIGNAL = SCALE * SIGNAL;

fftMATRIX = fft(eye(DATA_LENGTH))/sqrt(DATA_LENGTH);

DOUT = zeros(PATHS,DATA_LENGTH);
INTEGRATOR = zeros(PATHS,DATA_LENGTH);
RESIDUE = zeros(PATHS,DATA_LENGTH);


for ii = 1:DATA_LENGTH*WARMUP_CYCLES+1
    RESIDUE(:,1) = SIGNAL(mod(ii-1,DATA_LENGTH)+1)-PHI(:,mod(ii-1,DATA_LENGTH)+1).*DOUT(:,1);
    INTEGRATOR(:,1) = INTEGRATOR(:,1)+PHI(:,mod(ii-1,DATA_LENGTH)+1).*RESIDUE(:,1);
%    DOUT(:,1) = floor(INTEGRATOR(:,1));
    DOUT(:,1) = 2*(INTEGRATOR(:,1)>0)-1;
%     RESDIUE(1,1)
end

for ii = 2:DATA_LENGTH
    RESIDUE(:,ii) = SIGNAL(ii)-PHI(:,ii).*DOUT(:,ii-1);
    INTEGRATOR(:,ii) = PHI(:,ii).*RESIDUE(:,ii)+INTEGRATOR(:,ii-1);
%    DOUT(:,1) = floor(INTEGRATOR(:,1));
    DOUT(:,ii) = 2*(INTEGRATOR(:,ii)>0)-1;
end

Y = sum(DOUT,2);
A = PHI*(fftMATRIX^-1);

tic;
alg_choice = setOpts('alg_choice',1);
if alg_choice == 1
    BP_OPTIONS = spgSetParms(...
    'verbosity',    setOpts('bp_verbosity',1),...       % 0=quiet, 1=some output, 2=more output
    'iterations',   setOpts('bp_iterations',1e07),...    % Max number of iterations                              (default is 10*m)
    'nPrevVals',    setOpts('bp_nPrevVals',3),...      % Number previous func values for linesearch            (default is 3)
    'bpTol',        setOpts('bp_bpTol',1e-08),...   % Tolerance for basis pursuit solution                  (default is 1e-8)
    'optTol',       setOpts('bp_optTol',1e-06),...   % Optimality tolerance                                  (default is 1e-6)
    'decTol',       setOpts('bp_decTol',1e-06),...   % Req'd rel. change in primal obj. for Newton           (default is 1e-6)
    'stepMin',      setOpts('bp_stepMin',1e-16),...   % Minimum spectral step                                 (default is 1e-16)
    'stepMax',      setOpts('bp_stepMax',1e+05),...   % Maximum spectral step                                 (default is 1e5)
    'subspaceMin',  setOpts('bp_subspaceMin',0));        % 0=no subspace minimization, 1=subspace minimization.  (default is 0)
    [x,r] = spg_bp(A,Y,BP_OPTIONS);
elseif alg_choice == 2
    k = setOpts('OMP_residual',1e-6);
    [x,r] = OMP( A, Y, k);
elseif alg_choice == 3
    k = setOpts('CoSaMP_estimate_sparsity',5);
    [x,r] = CoSaMP( A, Y, k);
else
    error('alg_choice must be either 1 for sgp_bp, or 2 for OMP, or 3 for CoSaMP');
end
X_STAR = (x+randn(DATA_LENGTH,1)/10^15)/SCALE;
SIGNAL_STAR = real((fftMATRIX^-1)*(X_STAR));

figure(1),
plot(1:numPointsFFT,20*log10(abs(fftMATRIX*SIGNAL')));
xlim([0 numPointsFFT/2]);
title('Input Signal FFT Graph');
xlabel('Frequency');
ylabel('dB');

figure(2)
subplot(2,1,1), 
plot(...
    1:numPointsFFT,20*log10(abs(X_STAR)),'b',...
    1:numPointsFFT,20*log10(abs(floor(fftMATRIX*SIGNAL'))),'r');
axis([0,numPointsFFT/2,-5,max(20*log10(abs(floor(fftMATRIX*SIGNAL'))))+5]);
title('Compare Original Signal and Reconstructed Signal in Frequency Domain');
xlabel('Frequency');
ylabel('dB');

subplot(2,1,2), plot(1:numPointsFFT,20*log10(abs(abs(X_STAR)-abs(fftMATRIX*SIGNAL'))));
xlim([0 numPointsFFT/2]);
title('Difference between Original Signal and Reconstructed Signal in Frequency Domain');
xlabel('Frequency');
ylabel('dB');
axis([0,numPointsFFT,-20,max(20*log10(abs(abs(X_STAR)-abs(fftMATRIX*SIGNAL'))))+5]);


figure(3)
subplot(2,1,1), plot(1:length(SIGNAL),SIGNAL,1:length(SIGNAL),SIGNAL_STAR);
xlabel('Time');
ylabel('Amplitude');
xlim([0 numPointsFFT]);
title('Compare Original Signal and Reconstructed Signal FFT Graph in Time Domain');

subplot(2,1,2), plot(1:length(SIGNAL),SIGNAL'-SIGNAL_STAR);
xlabel('Time');
ylabel('Error in Amplitude');
xlim([0 numPointsFFT]);
title(sprintf('Difference in Time Domain, Total Residual in Power is %0.4e',sum(abs(r).^2)));

end