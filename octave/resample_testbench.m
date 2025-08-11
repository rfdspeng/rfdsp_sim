% Test bench for sanity-checking resample_custom.m

clear; clc; close all;

rand("seed","reset");
addpath("tools");
addpath("models");

en_plot = 1;
en_tprecode = 1;

p = 1+round(rand()*15); % 1 to 16
q = 1+round(rand()*15); % 1 to 16

disp("--- Resample Ratio: P / Q ---");
disp(p); disp(q); disp("");

% Waveform
nsym = 120; bw = 20; scs = 15; num_sc = 1200; start_sc = 600-num_sc/2; modorder = 4;
ncp = 10; wola = 2;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola);
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% Resample and inverse resample
[x1,cfg1] = resample_custom(x,p,q);
[x2,cfg2] = resample_custom(x1,q,p);

% Plot PSDs
if en_plot
  fs1 = fs*p/q;
  fsmax = max(fs,fs1);
  rbw = scs/1000/2^6;
  [p,f] = calculate_psd(x,fs,rbw); p = scale_psd(p,f,bw,scs,start_sc,num_sc);
  [p1,f1] = calculate_psd(x1,fs1,rbw); p1 = scale_psd(p1,f1,bw,scs,start_sc,num_sc);
  p2 = calculate_psd(x2,fs,rbw); p2 = scale_psd(p2,f,bw,scs,start_sc,num_sc);

  figure; hold on;
  plot(f,10*log10(p),'linewidth',5);
  plot(f1,10*log10(p1),'linewidth',3);
  plot(f,10*log10(p2),'linewidth',1);
  title("PSD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  legend("X","X1","X2",'location','south');
  xlim([-fsmax/2 fsmax/2]);
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x(1+wola_len/2:end);
y_evm = x2(1+wola_len/2:end);
cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);

% Output statements
disp("--- First Filter: n / f / a / w / Max Ripple (dB) / Min Atten (dB) ---");
disp(cfg1.n); disp(cfg1.f); disp(cfg1.a); disp(cfg1.w); disp(cfg1.max_ripple); disp(cfg1.min_atten); disp("");

disp("--- Second Filter: n / f / a / w / Max Ripple (dB) / Min Atten (dB) ---");
disp(cfg2.n); disp(cfg2.f); disp(cfg2.a); disp(cfg2.w); disp(cfg2.max_ripple); disp(cfg2.min_atten); disp("");

disp("--- EVM (%) / SNR (dB) After Resampling Twice ---");
disp(evm); disp(snr); disp("");
