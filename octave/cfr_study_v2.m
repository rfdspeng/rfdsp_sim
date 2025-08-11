% Test bench for experimenting with CFR techniques

clear; clc; close all;

addpath("models");
addpath("tools");

en_plot = 1;

en_tprecode = 0;
target_papr = 6;
osr = 4; % upsampling prior to CFR

% Waveform
bitwidth = 16; setpoint = -16;
nsym = 120; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-num_sc/2; modorder = 4;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% Upsampling
x = resample_custom(x,osr,1);
fs = fs*osr;

% Scaling
x = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1);
papr = calculate_papr(x,99.99);

% Hard clipping
th = 10^((setpoint+target_papr)/20)*2^(bitwidth-1);
env = abs(x);
ph = angle(x);
env(env > th) = th;
x_hc = env.*exp(1j*ph);
papr_hc = calculate_papr(x_hc,99.99);

% Plot hard clipping transfer curve
if en_plot
  figure; plot(abs(x),abs(x_hc),'.');
  title("Hard Clipping"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
endif

% Calculate hard clipping EVM
x_evm = resample_custom(x,1,osr);
x_evm = x_evm(1+wola_len/2:end);
y_evm = resample_custom(x_hc,1,osr);
y_evm = y_evm(1+wola_len/2:end);
cfg_evm.en_plot = en_plot; cfg_evm.title = "Hard Clipping";
evm_hc = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_hc = -20*log10(evm_hc/100);

% Calculate hard clipping ACLR
cfg_aclr = cfg_evm; cfg_aclr.fs = fs;
[aclr,aclrm,aclrp] = aclr_calculator(cfg_aclr,x);
[aclr_hc,aclrm_hc,aclrp_hc] = aclr_calculator(cfg_aclr,x_hc);

% Plot hard clipping PSDs
if en_plot
  rbw = scs/1000/2^2;
  [px,f] = calculate_psd(x,fs,rbw); px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(x_hc,fs,rbw); py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',4); plot(f,10*log10(py),'linewidth',1);
  title("Hard Clipping"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Output metrics
disp("--- PAPR: Base / HC ---");
disp(papr); disp(papr_hc); disp("");

disp("--- EVM: HC ---");
disp(evm_hc); disp("");

disp("--- ACLR: Base / HC ---");
disp(aclr); disp(aclr_hc); disp("");
