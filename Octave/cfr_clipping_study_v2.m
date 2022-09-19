% Test bench for clipping and filtering

clear; clc; close all;

addpath("models");
addpath("tools");

en_plot = 1;

en_tprecode = 0;
num_sc = 1200;
start_sc = 600-num_sc/2;
target_papr = 6;
niter = 10;

%{
en_tprecode = 1;
num_sc = 12;
start_sc = 0;
target_papr = 3.5;
%}

osr = 4; % upsampling prior to CFR

% Waveform
bitwidth = 16; setpoint = -16;
nsym = 120; bw = 20; scs = 15; modorder = 4;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% Upsampling
x_evm = x(1+wola_len/2:end);
x = resample_custom(x,osr,1);
fs = fs*osr;

% Scaling
x = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1);
papr = calculate_papr(x,99.99);
[ccdf,bins] = calculate_ccdf(x);
cfg_aclr = cfg_evm; cfg_aclr.fs = fs;
[aclr,aclrm,aclrp] = aclr_calculator(cfg_aclr,x);

% Filter coefficients
nrb = bw*5; obw = nrb*12*scs/1000; gb = (bw-obw)/2;
fpass = obw/fs; % normalized passband
fstop = (bw+gb)/fs; % normalized stopband
ripple = 0.01; atten = 5; % dB specs
n = 120; en_check = 1;
[b,cfgb] = firls_wrapper(n,fpass,fstop,ripple,atten,en_check);

% Iterative clipping and filtering
target_paprs = linspace(target_papr,papr,niter-1);
delta = min(diff(target_paprs));
target_paprs = [target_papr-delta target_paprs];
target_paprs = fliplr(target_paprs);
paprs = zeros(1,niter);
evms = zeros(1,niter);
snrs = zeros(1,niter);
aclrs = zeros(1,niter);
ccdfs = zeros(niter,length(ccdf));
x_hcs = zeros(niter,length(x));
x_hcfilts = zeros(niter,length(x));
xiter = x;
for idx = 1:niter
  target_papr_i = target_paprs(idx);

  % Clipping
  th = 10^((setpoint+target_papr_i)/20)*2^(bitwidth-1);
  env = abs(xiter);
  ph = angle(xiter);
  env(env > th) = th;
  x_hc = env.*exp(1j*ph);

  % Filtering
  x_hcfilt = filter(b,1,[x_hc zeros(1,n/2)]);
  x_hcfilt = x_hcfilt(1+n/2:end);

  % EVM
  y_evm = resample_custom(x_hcfilt,1,osr);
  y_evm = y_evm(1+wola_len/2:end);
  cfg_evm.en_plot = 0; %cfg_evm.title = "HC";
  evms(idx) = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
  snr = -20*log10(evms(idx)/100);

  % ACLR
  aclrs(idx) = aclr_calculator(cfg_aclr,x_hcfilt);

  % Prepare for next iteration
  xiter = x_hcfilt;
  x_hcs(idx,:) = x_hc;
  x_hcfilts(idx,:) = x_hcfilt;
  paprs(idx) = calculate_papr(x_hcfilt,99.99);
  ccdfs(idx,:) = calculate_ccdf(x_hcfilt);
endfor

% Plot CCDFs
if en_plot && 0
  figure; hold on;
  semilogy(bins,ccdf);
  for idx = 1:niter
    semilogy(bins,ccdfs(idx,:));
  endfor
endif

% Plot AMAM
if en_plot && 0
  figure; hold on;
  for idx = 1:niter
    plot(abs(x),abs(x_hcfilts(idx,:)),'.');
  endfor
  title("CFR AMAM"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  %legend("HC","HC+Filt",'location','south');
  grid on; set(gca,"fontsize",40);
endif

% Plot PSDs
if en_plot && 0
  rbw = scs/1000/2^2;
  [p,f] = calculate_psd(x,fs,rbw); p = scale_psd(p,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(p),'linewidth',5);
  for idx = 1:niter
    p_hcfilt = calculate_psd(x_hcfilts(idx,:),fs,rbw); p_hcfilt = scale_psd(p_hcfilt,f,bw,scs,start_sc,num_sc);
    plot(f,10*log10(p_hcfilt),'linewidth',1);
  endfor
  title("CFR PSDs"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  %legend("Mod","HC","HC+Filt",'location','south');
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",40);
endif

% Output metrics
disp("--- PAPR: Base / Iterations ---");
disp(papr);
for idx = 1:niter
  disp(paprs(idx));
endfor
disp("");

disp("--- EVM: Iterations ---");
for idx = 1:niter
  disp(evms(idx));
endfor
disp("");

disp("--- ACLR: Base / Iterations ---");
disp(aclr);
for idx = 1:niter
  disp(aclrs(idx));
endfor
disp("");
