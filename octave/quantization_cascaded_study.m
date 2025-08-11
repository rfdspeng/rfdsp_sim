% Simulate how quantization noise adds through a signal processing chain
% 1: Sanity-check SQNR after quantizing from floating point to fixed point
% 2: Run fixed-point waveform through fixed-point LPF multiple times and check SQNR

clear; clc; close all;

addpath("models");
addpath("tools");

en_fd_eq = 1;
en_plot = 1;
niter = 5;
rs = 10; %rs = 24;

% Generate floating-point and fixed-point waveforms
bitwidth = 16; setpoint = -16;
nsym = 120*5; bw = 20; scs = 15; num_sc = 100; start_sc = 600-num_sc/2; modorder = 4; en_tprecode = 0;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
x_fl = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1); % set waveform to correct setpoint and bitwidth
x = round(x_fl); % fixed-point conversion

% TD SNR after fixed-point conversion
evm_td_1 = sqrt((x_fl-x)*(x_fl-x)'/(x_fl*x_fl'))*100;
snr_td_1 = -20*log10(evm_td_1/100);

% FD SNR after fixed-point conversion
x_evm = x_fl(cfg_evm.wola_len/2+1:end);
y_evm = x(cfg_evm.wola_len/2+1:end);
cfg_evm.en_plot = en_plot; cfg_evm.en_fd_eq = en_fd_eq;
cfg_evm.title = "Floating Point to Fixed Point";
evm_1 = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_1 = -20*log10(evm_1/100);

% Expected FD and TD SNR after fixed-point conversion
snr_1_expected = 6.02*bitwidth+1.76+(setpoint+0)+10*log10(fs/(num_sc*scs/1000));
snr_td_1_expected = 6.02*bitwidth+1.76+(setpoint+0);

% Plot floating-point and fixed-point PSDs
if en_plot
  rbw = scs/1000/2^3;
  [p,f] = calculate_psd(x,fs,rbw);
  p_fl = calculate_psd(x_fl,fs,rbw);
  p = scale_psd (p,f,bw,scs,start_sc,num_sc);
  p_fl = scale_psd (p_fl,f,bw,scs,start_sc,num_sc);
  figure; hold on; plot(f,10*log10(p_fl),'linewidth',10); plot(f,10*log10(p),'linewidth',7);
  title("Floating Point to Fixed Point");
  xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  set(gca,"ytick",-120:3:5,"fontsize",20);
  legend("Floating Point","Fixed Point");
  xlim([-fs/2 fs/2]); ylim([-120 5]);
  grid on;
endif

% Generate LPF coefficients
n = 100;
fpass = 0.7; fstop = 0.8;
f_lpf = [0 fpass fstop 1];
a = [1 1 0 0];
dev = [1-10^(-0.01/20) 10^(-50/20)];
w = max(dev)./dev; w = w.*w;
b = firls(n,f_lpf,a,w); b = b/sum(b);
[h,om] = freqz(b,1,2^16); om = om/pi;
hpass = h(om <= fpass);
ripple_fl = max(abs(20*log10(abs(hpass))));
b_fl = b; h_fl = h;

% Quantize LPF coefficients
b = round(b*2^rs); %b = b*2^5;
[h,om] = freqz(b/sum(b),1,2^16); om = om/pi;
hpass = h(om <= fpass);
ripple = max(abs(20*log10(abs(hpass))));

% Plot floating-point and fixed-point LPF response
if en_plot
  figure; hold on;
  plot(om*fs/2,20*log10(abs(h_fl)),'linewidth',5);
  plot(om*fs/2,20*log10(abs(h)),'linewidth',5);
  xlim([0 fs/2]); ylim([-120 2]); grid on;
  title("LPF Magnitude Response");
  xlabel("Frequency (MHz)"); ylabel("Magnitude Response (dB)");
  set(gca,"fontsize",20);
  legend("Floating-Point Coefficients","Fixed-Point Coefficients");
endif

% Run waveform through LPF niter times
if en_plot, figure; hold on; plot(f,10*log10(p)); endif % plot fixed-point PSD after every filtering
x_iter = x; % fixed point: Quantize after every filter
x_iter_fl = x_fl; % floating-point reference
for idx = 1:niter
  x_iter_fl = filter(b/2^rs,1,[x_iter_fl zeros(1,n/2)]);
  x_iter_fl = x_iter_fl(1+n/2:end);

  x_iter = filter(b/2^rs,1,[x_iter zeros(1,n/2)]);
  x_iter = x_iter(1+n/2:end);
  x_iter = round(x_iter);

  if en_plot
    p_iter = calculate_psd(x_iter,fs,rbw);
    p_iter = scale_psd(p_iter,f,bw,scs,start_sc,num_sc);
    plot(f,10*log10(p_iter));
  endif
endfor
if en_plot
  title("Fixed-Point Cascade");
  xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  set(gca,"ytick",-120:3:5,"fontsize",20);
  xlim([-fs/2 fs/2]); ylim([-120 5]);
  grid on;
endif

% Compare PSDs
% 1 - Fixed-point input
% 2 - Fixed-point output
% 3 - Floating-point output
if en_plot
  figure; hold on;
  p_iter_fl = calculate_psd(x_iter_fl,fs,rbw);
  p_iter_fl = scale_psd(p_iter_fl,f,bw,scs,start_sc,num_sc);
  plot(f,10*log10(p),'linewidth',10); % 1
  plot(f,10*log10(p_iter),'linewidth',7); % 2
  plot(f,10*log10(p_iter_fl),'linewidth',5); % 3
  title("Input and Output After Iterations");
  set(gca,"ytick",-120:3:5,"fontsize",20);
  legend("Fixed-Point Input","Fixed-Point Output","Floating-Point Reference");
  xlim([-fs/2 fs/2]); ylim([-120 5]);
  grid on;
endif

%{
% Floating-point SNR after iterations
x_evm = x(cfg_evm.wola_len/2+1:end);
y_evm = x_iter_fl(cfg_evm.wola_len/2+1:end);
cfg_evm.en_plot = 0; cfg_evm.en_fd_eq = en_fd_eq;
cfg_evm.title = "Floating-Point Iterations";
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_iter_fl = -20*log10(evm/100);

% Fixed-point SNR after iterations
x_evm = x(cfg_evm.wola_len/2+1:end);
y_evm = x_iter(cfg_evm.wola_len/2+1:end);
cfg_evm.en_plot = 0; cfg_evm.en_fd_eq = en_fd_eq;
cfg_evm.title = "Fixed-Point Iterations";
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_iter = -20*log10(evm/100);
%}

% Fixed-point SNR after iterations, with LPF response decoupled
x_evm = x_iter_fl(cfg_evm.wola_len/2+1:end);
y_evm = x_iter(cfg_evm.wola_len/2+1:end);
cfg_evm.en_plot = en_plot; cfg_evm.en_fd_eq = en_fd_eq;
cfg_evm.title = "Fixed-Point Iterations";
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_iter = -20*log10(evm/100);

% TD SNR after iterations, with LPF response decoupled
evm = sqrt((x_iter_fl-x_iter)*(x_iter_fl-x_iter)'/(x_iter_fl*x_iter_fl'))*100;
snr_td_iter = -20*log10(evm/100);

% Final results
disp("--- SQNR: FD Calculated / Expected / TD Calculated / Expected ---");
disp(snr_1); disp(snr_1_expected); disp(snr_td_1); disp(snr_td_1_expected); disp("");

disp("--- LPF In-Band Ripple: Floating / Fixed ---");
disp(ripple_fl); disp(ripple); disp("");

disp("--- SQNR After Iterations: FD / TD ---");
disp(snr_iter); disp(snr_td_iter); disp("");

disp("--- Signal RMS: Input / Floating Output / Fixed Output ---");
disp(rms(x)); disp(rms(x_iter_fl)); disp(rms(x_iter)); disp("");


%snr_td_expected = 6.02*bitwidth+1.76+(setpoint+0)+10*log10(2/fpass/2);
%disp("---Expected TD SNR (dB) after floating- to fixed-point quantization---");

% The quantization noise floor growth is affected by the coefficient bitwidth.
% If I keep the full floating-point resolution for the coefficients, even if I round after the filtering,
% the noise floor does not grow.

% Is this because some coefficients are too small? So the overall gain for some taps is very low, so you get more quantization error
% Quantizing from B1 to B is one thing, but quantizing after multiplying with coefficients seems different
% If you increase the filter gain, then results are good

% Somehow FD SNR actually improves if I increase rs
