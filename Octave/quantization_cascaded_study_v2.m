% Simulate how quantization noise adds through a signal processing chain

clear; clc; close all;

en_plot = 1;

bitwidth = 16; setpoint = -16; % data params

nsym = 120*5; bw = 20; scs = 15; num_sc = 100; start_sc = 600-num_sc/2; modorder = 4; en_tprecode = 0;
ncp = 10; wola = 2; seed = 1;

[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
x_fl = x/rms(x)*10^(setpoint/20)*2^(bitwidth-1);
x = round(x_fl);
evm_td = sqrt((x_fl-x)*(x_fl-x)'/(x_fl*x_fl'))*100;
snr_td = -20*log10(evm_td/100);

x_evm = x_fl(cfg_evm.wola_len/2+1:end);
y_evm = x(cfg_evm.wola_len/2+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
%disp(evm);
disp("---FD SNR (dB) after floating- to fixed-point quantization---");
disp(snr);
%disp(evm_td);
disp("---TD SNR (dB) after floating- to fixed-point quantization---");
disp(snr_td);

snr_expected = 6.02*bitwidth+1.76+(setpoint+0)+10*log10(fs/(num_sc*scs/1000));
snr_td_expected = 6.02*bitwidth+1.76+(setpoint+0);
disp("---Expected FD SNR (dB) after floating- to fixed-point quantization---");
disp(snr_expected);
disp("---Expected TD SNR (dB) after floating- to fixed-point quantization---");
disp(snr_td_expected);

rbw = scs/1000/2^3;
[p,f] = calculate_psd(x,fs,rbw);
p_fl = calculate_psd(x_fl,fs,rbw);
p = scale_psd (p,f,bw,scs,start_sc,num_sc);
p_fl = scale_psd (p_fl,f,bw,scs,start_sc,num_sc);
figure; hold on; plot(f,10*log10(p)); plot(f,10*log10(p_fl));
set(gca,"ytick",-120:3:5,"fontsize",20);
xlim([-fs/2 fs/2]); ylim([-120 5]);
grid on;

% Generate low-pass filter coefficients
n = 100;
fpass = 0.7; fstop = 0.8;
f = [0 fpass fstop 1];
a = [1 1 0 0];
dev = [1-10^(-0.01/20) 10^(-50/20)];
w = max(dev)./dev; w = w.*w;
b = firls(n,f,a,w); b = b/sum(b);
[h,om] = freqz(b,1,2^16); om = om/pi;
hpass = h(om <= fpass);
ripple = max(abs(20*log10(abs(hpass))));
disp("---Passband ripple (dB)---");
disp(ripple);
b_fl = b; h_fl = h;

% Quantize coefficients
rs = 10;
b = round(b*2^rs);
[h,om] = freqz(b/sum(b),1,2^16); om = om/pi;
hpass = h(om <= fpass);
ripple = max(abs(20*log10(abs(hpass))));
disp("---Passband ripple (dB) after coefficients are quantized---");
disp(ripple);

figure; hold on;
plot(om*fs/2,20*log10(abs(h_fl)),'linewidth',5);
plot(om*fs/2,20*log10(abs(h)),'linewidth',5);
xlim([0 fs/2]); ylim([-120 2]); grid on;
set(gca,"fontsize",20);

% Run waveform through low-pass filter multiple times
niter = 9;
figure; hold on;
rbw = scs/1000/2^3;
[p,f] = calculate_psd(x,fs,rbw);
p_fl = calculate_psd(x_fl,fs,rbw);
p = scale_psd(p,f,bw,scs,start_sc,num_sc);
p_fl = scale_psd(p_fl,f,bw,scs,start_sc,num_sc);
plot(f,10*log10(p_fl)); plot(f,10*log10(p));
x_iter = x;
x_iter_fl = x; % no quantization after filtering (input data is still quantized)
for idx = 1:niter
  %x_iter = filter(b,1,[x_iter zeros(1,n/2)]);
  %x_iter = round(x_iter/2^rs);

  x_iter = filter(b/2^rs,1,[x_iter zeros(1,n/2)]);
  x_iter = x_iter(1+n/2:end);
  x_iter = round(x_iter);
  %_iter = floor(x_iter);

  x_iter_fl = filter(b/2^rs,1,[x_iter_fl zeros(1,n/2)]);
  x_iter_fl = x_iter_fl(1+n/2:end);

  p_iter = calculate_psd(x_iter,fs,rbw);
  p_iter = scale_psd(p_iter,f,bw,scs,start_sc,num_sc);
  plot(f,10*log10(p_iter));
endfor
set(gca,"ytick",-120:3:5,"fontsize",20);
xlim([-fs/2 fs/2]); ylim([-120 5]);
grid on;

% Compare floating point to fixed point
figure; hold on;
p_iter_fl = calculate_psd(x_iter_fl,fs,rbw);
p_iter_fl = scale_psd(p_iter_fl,f,bw,scs,start_sc,num_sc);
plot(f,10*log10(p)); % quantized input
plot(f,10*log10(p_iter)); % quantized output
plot(f,10*log10(p_iter_fl)); % un-quantized output
set(gca,"ytick",-120:3:5,"fontsize",20);
xlim([-fs/2 fs/2]); ylim([-120 5]);
grid on;

% Calculate EVM
x_evm = x(cfg_evm.wola_len/2+1:end);
y_evm = x_iter_fl(cfg_evm.wola_len/2+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
disp("---FD SNR (dB) after running through LPF without output quantization---");
disp(snr);

x_evm = x(cfg_evm.wola_len/2+1:end);
y_evm = x_iter(cfg_evm.wola_len/2+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
disp("---FD SNR (dB) after running through LPF with output quantization---");
disp(snr);

x_evm = x_iter_fl(cfg_evm.wola_len/2+1:end);
y_evm = x_iter(cfg_evm.wola_len/2+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
disp("---FD SNR (dB) after running through LPF with output quantization with filter response decoupled---");
disp(snr);

evm_td = sqrt((x_iter_fl-x_iter)*(x_iter_fl-x_iter)'/(x_iter_fl*x_iter_fl'))*100;
snr_td = -20*log10(evm_td/100);
disp(snr_td);
%snr_td_expected = 6.02*bitwidth+1.76+(setpoint+0)+10*log10(2/fpass/2);
%disp("---Expected TD SNR (dB) after floating- to fixed-point quantization---");

% The quantization noise floor growth is affected by the coefficient bitwidth.
% If I keep the full floating-point resolution for the coefficients, even if I round after the filtering,
% the noise floor does not grow.
