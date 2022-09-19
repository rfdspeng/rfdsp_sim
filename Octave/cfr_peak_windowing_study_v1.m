% Test bench for peak windowing

clear; clc; close all;

addpath("models");
addpath("tools");

en_plot = 1;

en_tprecode = 0;
num_sc = 1200;
start_sc = 600-num_sc/2;
target_papr = 6;

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
%x = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1);
x = x/max(abs(x));
setpoint = 20*log10(rms(x));

% Peak windowing
%th = 10^((setpoint+target_papr)/20)*2^(bitwidth-1);
th = 10^((setpoint+target_papr)/20);
env = abs(x); ph = angle(x);
pdx = env > th;
s_hc = ones(size(x)); % clipping coefficients
s_hc(pdx) = th./env(pdx);
%figure; hold on; plot(env); plot(s_hc);


winlen = 15;
win = hanning(winlen);
win = ones(1,winlen)/winlen*(winlen-1)/1.2;
%win1 = upsample(win(2:end-1),2);
%win = win1(1:end-1); win = win.';
win1 = upsample(win(2:end-1),3);
win = win1(1:end-2); win = win.';
winlen = length(win);
d = (winlen-1)/2;
%}

s_pw = filter(win,1,[1-s_hc zeros(1,d)]);
s_pw = s_pw(1+d:end);
s_pw = 1-s_pw;
disp(sum(s_pw > s_hc));
disp(sum(s_pw < 0));
%figure; hold on; plot(s_hc); plot(s_pw);

s_pw(s_pw < 0) = 0;


% Window
x_hc = x.*s_hc;
x_pw = x.*s_pw;

% PAPR
papr = calculate_papr(x,99.99);
papr_hc = calculate_papr(x_hc,99.99);
papr_pw = calculate_papr(x_pw,99.99);

% CCDF
[ccdf,bins] = calculate_ccdf(x);
ccdfhc = calculate_ccdf(x_hc);
ccdfpw = calculate_ccdf(x_pw);

if en_plot && 0
  figure; hold on;
  semilogy(bins,ccdf);
  semilogy(bins,ccdfhc);
  semilogy(bins,ccdfpw);

endif

% AMAM
if en_plot && 0
  figure; hold on;
  plot(abs(x),abs(x_hc),'.');
  plot(abs(x),abs(x_pw),'.');
  title("CFR AMAM"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  legend("HC","PW",'location','south');
  grid on; set(gca,"fontsize",40);
endif

% EVM
x_evm = resample_custom(x,1,osr);
x_evm = x_evm(1+wola_len/2:end);
y_evm = resample_custom(x_hc,1,osr);
y_evm = y_evm(1+wola_len/2:end);
cfg_evm.en_plot = en_plot && 0; cfg_evm.title = "HC";
evm_hc = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_hc = -20*log10(evm_hc/100);

y_evm = resample_custom(x_pw,1,osr);
y_evm = y_evm(1+wola_len/2:end);
cfg_evm.title = "PW";
evm_pw = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_pw = -20*log10(evm_pw/100);

% ACLR
cfg_aclr = cfg_evm; cfg_aclr.fs = fs;
[aclr,aclrm,aclrp] = aclr_calculator(cfg_aclr,x);
[aclr_hc,aclrm_hc,aclrp_hc] = aclr_calculator(cfg_aclr,x_hc);
[aclr_pw,aclrm_pw,aclrp_pw] = aclr_calculator(cfg_aclr,x_pw);

% Signal PSDs
if en_plot
  rbw = scs/1000/2^2;
  [p,f] = calculate_psd(x,fs,rbw); p = scale_psd(p,f,bw,scs,start_sc,num_sc);
  p_hc = calculate_psd(x_hc,fs,rbw); p_hc = scale_psd(p_hc,f,bw,scs,start_sc,num_sc);
  p_pw = calculate_psd(x_pw,fs,rbw); p_pw = scale_psd(p_pw,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(p),'linewidth',5);
  plot(f,10*log10(p_hc),'linewidth',3);
  plot(f,10*log10(p_pw),'linewidth',1);
  title("CFR PSDs"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  legend("Mod","HC","PW",'location','south');
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",40);
endif

% Window PSDs
if en_plot
  rbw = scs/1000/2^2;
  [w_hc,f] = calculate_psd(s_hc,fs,rbw); %w_hc = scale_psd(w_hc,f,bw,scs,start_sc,num_sc);
  w_pw = calculate_psd(s_pw,fs,rbw); %w_pw = scale_psd(w_pw,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(w_hc),'linewidth',3);
  plot(f,10*log10(w_pw),'linewidth',1);
  title("Window PSDs"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  legend("HC","PW",'location','south');
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",40);
endif

% Output metrics
disp("--- PAPR: Base / HC / PW ---");
disp(papr); disp(papr_hc); disp(papr_pw); disp("");

disp("--- EVM: HC / PW ---");
disp(evm_hc); disp(evm_pw); disp("");

disp("--- ACLR: Base / HC / PW ---");
disp(aclr); disp(aclr_hc); disp(aclr_pw); disp("");
