% Derive RF amplifier model from IIP2, IIP3
% Verify the model using two-tone test
% Run OFDM waveform through model

clear; clc; close all;

addpath("tools");
addpath("models");

en_plot = 1;
niter = 1;

xlen = 2^17; nfft = round(xlen/2^6);

fbb = pi/32;
frf = pi/4;
%frf = pi/2; % the results are wonky if frf = pi/2 because some spurs wrap around and alias onto the spurs of interest

piip3 = rand()*50; % 0 to 50dBm
piip2 = piip3+rand()*30; % 0 to 80dBm
pin = -30-20*rand() + min(piip3,piip2); % -30 to -50dB lower than min(IIP3,IIP2)
pin_ofdm = pin + rand()*5-2.5; % OFDM signal power = tone power +/- 2.5dB
a1 = 1;

% In this simulation, pin is equal to the complex baseband signal power
% If pin is instead equal to the RF power, then you must add +3dB before converting to linear
r = 10^(pin/20);
r_ofdm = 10^(pin_ofdm/20);
a2 = a1/10^(piip2/20); %a2 = 0;
a3 = -4/3*a1*10^(-piip3/10);

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

% Post-processing
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
