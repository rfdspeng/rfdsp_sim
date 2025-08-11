% Verify RF and BB models for Tx LO FI RSB

clear; clc; close all;

addpath("models");
addpath("tools");

en_plot = 1;
en_fd_eq = 1;

gi = 1.1; gq = 0.9; phii = 2; phiq = 4;
%gi = 1; gq = 1; phii = 0; phiq = 0;

% Waveform
nsym = 120; bw = 20; scs = 15; modorder = 4; en_tprecode = 0;
ncp = 10; wola = 2; seed = 1;

num_sc = 50; start_sc = 600-num_sc-1; % RSB does not overlap signal
num_sc = 100; start_sc = 600-num_sc/2; % RSB overlaps signal
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% RF RSB model
fc = fs/4; wc = fc*2*pi/fs;
cfg_rsb.gi = gi; cfg_rsb.gq = gq;
cfg_rsb.phii = phii; cfg_rsb.phiq = phiq;
cfg_rsb.fc = fc; cfg_rsb.fs = fs;
y = tx_lo_fi_rsb_model_rf(cfg_rsb,x);

% Downconvert
len = length(y);
ybb = y.*exp(-1j*wc*(0:len-1));

% Filter 2wc
fpass = 0.3;
n = 50; f = [0 fpass 1-fpass 1];
a = [1 1 0 0];
b = firls(n,f,a);
[h,w] = freqz(b,1,2^17); w = w/pi;
nbw = num_sc*scs/1000/fs;
hpass = h(w <= nbw);
ripple = max(abs(20*log10(abs(hpass))));
disp(ripple);
ybb = filter(b,1,[ybb zeros(1,n/2)]);
ybb = ybb(1+n/2:end);

% Plot PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); %px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(y,fs,rbw); %py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  pybb = calculate_psd(ybb,fs,rbw);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',3); plot(f,10*log10(py),'linewidth',2); plot(f,10*log10(pybb),'linewidth',1);
  title("Tx LO FI RSB PSD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x(1+wola_len/2:end);
y_evm = ybb(1+wola_len/2:end);
cfg_evm.en_plot = en_plot;
cfg_evm.en_fd_eq = en_fd_eq;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
disp(evm);
disp(snr);

% BB RSB model
ybb = tx_lo_fi_rsb_model_bb(cfg_rsb,x);

% Plot PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); %px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  pybb = calculate_psd(ybb,fs,rbw);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',2); plot(f,10*log10(pybb),'linewidth',1);
  title("Tx LO FI RSB PSD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x(1+wola_len/2:end);
y_evm = ybb(1+wola_len/2:end);
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
disp(evm);
disp(snr);
