% Test bench for experimenting with CFR techniques

clear; clc; close all;

en_plot = 1;

target_papr = 5;
pw_len = 13;

% Waveform
bitwidth = 16; setpoint = -16;
nsym = 120; bw = 20; scs = 15; num_sc = 600; start_sc = 600-num_sc/2; modorder = 4; en_tprecode = 0;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
x = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1);
papr = calculate_papr(x,99.99);
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% Hard clipping
clip_th = 10^((setpoint+target_papr)/20)*2^(bitwidth-1);
peaks = abs(x) > clip_th;
c_hc = ones(size(x));
c_hc(peaks) = clip_th./abs(x(peaks));
x_hc = x.*c_hc;
papr_hc = calculate_papr(x_hc,99.99);

% Plot hard clipping transfer curve
if en_plot
  figure; plot(abs(x),abs(x_hc),'.');
  title("Hard Clipping"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
endif

% Calculate hard clipping EVM
x_evm = x(1+wola_len/2:end);
y_evm = x_hc(1+wola_len/2:end);
cfg_evm.en_plot = en_plot; cfg_evm.title = "Hard Clipping";
evm_hc = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_hc = -20*log10(evm_hc/100);

% Plot hard clipping PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(x_hc,fs,rbw); py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',4); plot(f,10*log10(py),'linewidth',1);
  title("Hard Clipping"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Peak windowing
win = hanning(pw_len).'; win = win/sum(win);
gd = (pw_len-1)/2;
c_pw = filter(win,1,[(1-c_hc) zeros(1,gd)]);
c_pw = c_pw(1+gd:end);
c_pw = 1-c_pw;
figure; hold on; plot(c_hc,'linewidth',10); plot(c_pw,'linewidth',7);
x_pw = x.*c_pw;
papr_pw = calculate_papr(x_pw,99.99);

% Plot peak windowing transfer curve
if en_plot
  figure; plot(abs(x),abs(x_pw),'.');
  title("Peak Windowing"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
endif

% Calculate peak windowing EVM
x_evm = x(1+wola_len/2:end);
y_evm = x_pw(1+wola_len/2:end);
cfg_evm.en_plot = en_plot; cfg_evm.title = "Peak Windowing";
evm_pw = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_pw = -20*log10(evm_pw/100);

% Plot peak windowing PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(x_pw,fs,rbw); py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',4); plot(f,10*log10(py),'linewidth',1);
  title("Peak Windowing"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Output metrics
disp("--- PAPR: Base / HC / PW ---");
disp(papr); disp(papr_hc); disp(papr_pw); disp("");

disp("--- EVM: HC / PW ---");
disp(evm_hc); disp(evm_pw); disp("");
