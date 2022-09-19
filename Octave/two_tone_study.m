% Derive RF amplifier model from IIP2, IIP3
% Verify the model using two-tone test
% Run OFDM waveform through model

clear; clc; close all;

addpath("tools");
addpath("models");

en_plot = 1;
niter = 1;

xlen = 2^17; nfft = round(xlen/2^6);
en_tprecode = 0;

fbb = pi/32;
frf = pi/4;
%frf = pi/2; % the results are wonky if frf = pi/2 because some spurs wrap around and alias onto the spurs of interest

rand("seed","reset");
%rand("state","reset");
%v = rand("state");
%disp("--- v ---");
%disp(v(1:3)); disp("");
piip3 = rand()*50; % 0 to 50dBm
piip2 = piip3+rand()*30; % 0 to 80dBm
piip2 = piip3+20+rand()*30;
pin = -30-20*rand() + min(piip3,piip2); % -30 to -50dB lower than min(IIP3,IIP2)
pin = -20-10*rand() + min(piip3,piip2);
pin_ofdm = pin + rand()*5-2.5; % OFDM signal power = tone power +/- 2.5dB
pin_ofdm = pin;
a1 = 1;

r = sqrt(2)*10^(pin/20);
r_ofdm = sqrt(2)*10^(pin_ofdm/20);
aiip2 = sqrt(2)*10^(piip2/20);
aiip3 = sqrt(2)*10^(piip3/20);
a2 = a1/aiip2;
a3 = -4/3*a1/aiip3^2;

if niter > 1, en_plot = 0; endif
im2s_bb = zeros(1,niter);
im2s_rf = zeros(1,niter);
im3s_bb = zeros(1,niter);
im3s_rf = zeros(1,niter);
for iterdx = 1:niter
  disp(iterdx);
  % Baseband input signal
  xbb = r*exp(1j*(fbb*(1:xlen)+rand()*2*pi)) + r*exp(1j*(-fbb*(1:xlen)+rand()*2*pi));

  % Baseband output signal
  ybb = a1*xbb + a2*abs(xbb).^2 + 3/4*a3*(abs(xbb).^2).*xbb;
  [ybb_psd,fbb_psd] = pwelch(ybb,hanning(nfft),[],nfft,2*pi,'centerdc',[],'none');

  % RF input signal
  lo = exp(1j*(frf*(1:xlen)+rand()*2*pi));
  xrf = real(xbb.*lo);

  % RF output signal
  yrf = a1*xrf + a2*xrf.^2 + a3*xrf.^3;
  [yrf_psd,frf_psd] = pwelch(yrf,hanning(nfft),[],nfft,2*pi,'centerdc',[],'none');

  % IM2 calc @ BB
  sl = -fbb-fbb/4; sh = -fbb+fbb/4;
  nl = -2*fbb-fbb/4; nh = -2*fbb+fbb/4;
  im2_bb_dbc = calculate_noise_dbc(ybb_psd,fbb_psd,sl,sh,nl,nh,en_plot,"IM2 @ BB");
  im2s_bb(iterdx) = im2_bb_dbc;

  % IM3 calc @ BB
  sl = -fbb-fbb/4; sh = -fbb+fbb/4;
  nl = -3*fbb-fbb/4; nh = -3*fbb+fbb/4;
  im3_bb_dbc = calculate_noise_dbc(ybb_psd,fbb_psd,sl,sh,nl,nh,en_plot,"IM3 @ BB");
  im3s_bb(iterdx) = im3_bb_dbc;

  % IM2 calc @ RF
  sl = +frf-fbb-fbb/4; sh = +frf-fbb+fbb/4;
  nl = -2*fbb-fbb/4; nh = -2*fbb+fbb/4;
  im2_rf_dbc = calculate_noise_dbc(yrf_psd,frf_psd,sl,sh,nl,nh,en_plot,"IM2 @ RF");
  im2s_rf(iterdx) = im2_rf_dbc;

  % IM3 calc @ RF
  sl = +frf-fbb-fbb/4; sh = +frf-fbb+fbb/4;
  nl = +frf-3*fbb-fbb/4; nh = +frf-3*fbb+fbb/4;
  im3_rf_dbc = calculate_noise_dbc(yrf_psd,frf_psd,sl,sh,nl,nh,en_plot,"IM3 @ RF");
  im3s_rf(iterdx) = im3_rf_dbc;
endfor

% Post-processing the two-tone simulation
sig_freqs_bb = [-fbb fbb];
spur_freqs_bb = [-3*fbb -2*fbb 2*fbb 3*fbb];
disp("--- BB Frequencies: Signal / Spurs ---");
disp(sig_freqs_bb); disp(spur_freqs_bb); disp("");

sig_freqs_rf = [-frf-fbb -frf+fbb frf-fbb frf+fbb];
spur_freqs_rf = [-3*frf-3*fbb -3*frf-fbb -3*frf+fbb -3*frf+3*fbb -frf-3*fbb -frf+3*fbb -2*fbb];
spur_freqs_rf = [spur_freqs_rf fliplr(abs(spur_freqs_rf))];
pdx = spur_freqs_rf > pi; ndx = spur_freqs_rf < -pi;
spur_freqs_rf(pdx) = spur_freqs_rf(pdx)-2*pi;
spur_freqs_rf(ndx) = spur_freqs_rf(ndx)+2*pi;
spur_freqs_rf = sort(spur_freqs_rf);
disp("--- RF Frequencies: Signal / Spurs ---");
disp(sig_freqs_rf); disp(spur_freqs_rf); disp("");

disp("--- Pin / Piip2 / Piip3 ---");
disp(pin); disp(piip2); disp(piip3); disp("");

im3_expected = 2*pin-2*piip3;
im3_bb_avg = mean(im3s_bb);
im3_rf_avg = mean(im3s_rf);
disp("--- IM3 (dBc): Expected / BB Sim / RF Sim ---");
disp(im3_expected); disp(im3_bb_avg); disp(im3_rf_avg); disp("");

im2_expected = pin-piip2;
im2_bb_avg = mean(im2s_bb);
im2_rf_avg = mean(im2s_rf);
disp("--- IM2 (dBc): Expected / BB Sim / RF Sim ---");
disp(im2_expected); disp(im2_bb_avg); disp(im2_rf_avg); disp("");

% OFDM waveform
nsym = 120; bw = 20; scs = 15; num_sc = 600; start_sc = 600-num_sc/2; modorder = 4;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
x = x/rms(x)*r_ofdm;
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% Nonlinear model
y = a1*x + a2*abs(x).^2 + 3/4*a3*(abs(x).^2).*x;

% Plot AMAM, AMPM
if en_plot
  figure; plot(abs(x),abs(y),'.');
  title("AMAM"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
  figure; plot(abs(x),180/pi*unwrap(angle(y)-angle(x)),'.');
  title("AMPM"); xlabel("Input Amplitude"); ylabel("Phase Modulation (Deg)");
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x(1+wola_len/2:end);
y_evm = y(1+wola_len/2:end);
cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);
disp("--- OFDM FD SNR (dB) ---");
disp(snr); disp("");

% Plot PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(y,fs,rbw); py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',4); plot(f,10*log10(py),'linewidth',1);
  title("PSD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif
