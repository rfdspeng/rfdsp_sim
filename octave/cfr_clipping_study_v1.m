% Test bench for clipping and filtering

clear; clc; close all;

addpath("models");
addpath("tools");

en_plot = 0;

en_tprecode = 0;
num_sc = 1200;
start_sc = 600-num_sc/2;
target_papr = 4.5;

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
x = resample_custom(x,osr,1);
fs = fs*osr;

% Scaling
x = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1);

% Clipping
th = 10^((setpoint+target_papr)/20)*2^(bitwidth-1);
env = abs(x);
ph = angle(x);
env(env > th) = th;
x_hc = env.*exp(1j*ph);

% Filtering
nrb = bw*5; obw = nrb*12*scs/1000; gb = (bw-obw)/2;
fpass = obw/fs; % normalized passband
fstop = (bw+gb)/fs; % normalized stopband
ripple = 0.01; atten = 5; % dB specs
n = 120; en_check = 1;
[b,cfgb] = firls_wrapper(n,fpass,fstop,ripple,atten,en_check);
x_hcfilt = filter(b,1,[x_hc zeros(1,n/2)]);
x_hcfilt = x_hcfilt(1+n/2:end);

% PAPR
papr = calculate_papr(x,99.99);
papr_hc = calculate_papr(x_hc,99.99);
papr_hcfilt = calculate_papr(x_hcfilt,99.99);

% CCDF
[ccdf,bins] = calculate_ccdf(x);
ccdfhc = calculate_ccdf(x_hc);
ccdfhcfilt = calculate_ccdf(x_hcfilt);

if en_plot
  figure; hold on;
  semilogy(bins,ccdf);
  semilogy(bins,ccdfhc);
  semilogy(bins,ccdfhcfilt);

endif

% AMAM
if en_plot
  figure; hold on;
  plot(abs(x),abs(x_hc),'.');
  plot(abs(x),abs(x_hcfilt),'.');
  title("CFR AMAM"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  legend("HC","HC+Filt",'location','south');
  grid on; set(gca,"fontsize",40);
endif

% EVM
x_evm = resample_custom(x,1,osr);
x_evm = x_evm(1+wola_len/2:end);
y_evm = resample_custom(x_hc,1,osr);
y_evm = y_evm(1+wola_len/2:end);
cfg_evm.en_plot = en_plot; cfg_evm.title = "HC";
evm_hc = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_hc = -20*log10(evm_hc/100);

y_evm = resample_custom(x_hcfilt,1,osr);
y_evm = y_evm(1+wola_len/2:end);
cfg_evm.title = "HC+Filt";
evm_hcfilt = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_hcfilt = -20*log10(evm_hcfilt/100);

% ACLR
cfg_aclr = cfg_evm; cfg_aclr.fs = fs;
[aclr,aclrm,aclrp] = aclr_calculator(cfg_aclr,x);
[aclr_hc,aclrm_hc,aclrp_hc] = aclr_calculator(cfg_aclr,x_hc);
[aclr_hcfilt,aclrm_hcfilt,aclrp_hcfilt] = aclr_calculator(cfg_aclr,x_hcfilt);

% PSDs
if en_plot
  rbw = scs/1000/2^2;
  [p,f] = calculate_psd(x,fs,rbw); p = scale_psd(p,f,bw,scs,start_sc,num_sc);
  p_hc = calculate_psd(x_hc,fs,rbw); p_hc = scale_psd(p_hc,f,bw,scs,start_sc,num_sc);
  p_hcfilt = calculate_psd(x_hcfilt,fs,rbw); p_hcfilt = scale_psd(p_hcfilt,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(p),'linewidth',5);
  plot(f,10*log10(p_hc),'linewidth',3);
  plot(f,10*log10(p_hcfilt),'linewidth',1);
  title("CFR PSDs"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  legend("Mod","HC","HC+Filt",'location','south');
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",40);
endif

% Output metrics
disp("--- PAPR: Base / HC / HC+Filt ---");
disp(papr); disp(papr_hc); disp(papr_hcfilt); disp("");

disp("--- EVM: HC / HC+Filt ---");
disp(evm_hc); disp(evm_hcfilt); disp("");

disp("--- ACLR: Base / HC / HC+Filt ---");
disp(aclr); disp(aclr_hc); disp(aclr_hcfilt); disp("");
