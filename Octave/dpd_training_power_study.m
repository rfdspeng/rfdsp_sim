% Test bench for experimenting with DPD training (EPT mode)

% Although DPD does lower the average output power compared to without DPD, it's still way better than simply reducing the PA input power -
% for the same output power reduction, DPD nearly completely recovers signal quality, but reducing the PA input power only improves EVM by
% 3dB.

clear; clc; close all;

addpath("tools");
addpath("models");

en_plot = 0;

% Waveform
nsym = 120; bw = 20; scs = 15; num_sc = 600; start_sc = 600-num_sc/2; modorder = 4; en_tprecode = 0;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
x = x/max(abs(x));
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% PA model
y = pa_model_v1(x);
comp = calculate_compression(x,y);

% Plot AMAM, AMPM
if en_plot
  figure; plot(abs(x),abs(y),'.');
  title("PA AMAM"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
  figure; plot(abs(x),180/pi*unwrap(angle(y)-angle(x)),'.');
  title("PA AMPM"); xlabel("Input Amplitude"); ylabel("Phase Modulation (Deg)");
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x(1+wola_len/2:end);
y_evm = y(1+wola_len/2:end);
cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr = -20*log10(evm/100);

% Plot PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(y,fs,rbw); py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',4); plot(f,10*log10(py),'linewidth',1);
  title("PA PSD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% DPD training
cfg_dpd.kernel_mask = [0 0 255];
                       %1 1 31;
                       %2 2 31
                       %3 3 31];
[c,cinfo] = dpd_training(cfg_dpd,x,y/max(abs(y)));
%[c,cinfo] = dpd_training(cfg_dpd,x,y);

% Predistort signal
env = abs(x);
x_dpd = 0;
for cdx = 1:length(c)
  delay = cinfo(cdx,1);
  env_delay = cinfo(cdx,2);
  env_order = cinfo(cdx,3);
  if delay >= 0, xl = [zeros(1,delay) x(1:end-delay)];
  else, xl = [x(1+delay:end) zeros(1,delay)]; endif
  if env_delay >= 0, envl = [zeros(1,env_delay) env(1:end-env_delay)];
  else, envl = [env(1+env_delay:end) zeros(1,env_delay)]; endif
  x_dpd = x_dpd + c(cdx)*xl.*envl.^env_order;
endfor

% Run predistorted signal through PA model
y_dpd = pa_model(x_dpd);

% Plot AMAM, AMPM
if en_plot
  figure; plot(abs(x),abs(y_dpd),'.');
  title("PA AMAM With DPD"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
  figure; plot(abs(x),180/pi*unwrap(angle(y_dpd)-angle(x)),'.');
  title("PA AMPM With DPD"); xlabel("Input Amplitude"); ylabel("Phase Modulation (Deg)");
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x(1+wola_len/2:end);
y_evm = y_dpd(1+wola_len/2:end);
cfg_evm.en_plot = en_plot;
evm_dpd = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_dpd = -20*log10(evm_dpd/100);

% Plot PSDs
if en_plot
  py_dpd = calculate_psd(y_dpd,fs,rbw); py_dpd = scale_psd(py_dpd,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',5); plot(f,10*log10(py),'linewidth',3); plot(f,10*log10(py_dpd),'linewidth',1);
  title("PA PSD With DPD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Calculate PAPR
[papr_x,avg_x,peak_x] = calculate_papr(x,100);
[papr_y,avg_y,peak_y] = calculate_papr(y,100);
[papr_x_dpd,avg_x_dpd,peak_x_dpd] = calculate_papr(x_dpd,100);
[papr_y_dpd,avg_y_dpd,peak_y_dpd] = calculate_papr(y_dpd,100);

%{
% How much does DPD reduce average power? (Should be the same as the compression)
disp("Input power without DPD minus input power with DPD");
disp(10*log10(avg_x/avg_x_dpd));
disp("Input peak power without DPD minus input peak power with DPD");
disp(10*log10(peak_x/peak_x_dpd));
disp("PA power without DPD minus PA power with DPD");
pin_backoff = avg_y/avg_y_dpd; % linear power backoff (for next experiment)
disp(10*log10(pin_backoff));
disp("PA peak power without DPD minus PA peak power with DPD");
disp(10*log10(peak_y/peak_y_dpd));
%}

% What happens if you just reduce input power?
pin_backoff = avg_y/avg_y_dpd;
pin_backoff = pin_backoff*3;
x_bo = x/sqrt(pin_backoff);

% Run backed off signal through PA model
y_bo = pa_model(x_bo);

% Plot AMAM, AMPM
if en_plot
  figure; plot(abs(x_bo),abs(y_bo),'.');
  title("PA AMAM With Backoff"); xlabel("Input Amplitude"); ylabel("Output Amplitude");
  grid on; set(gca,"fontsize",20);
  figure; plot(abs(x_bo),180/pi*unwrap(angle(y_bo)-angle(x_bo)),'.');
  title("PA AMPM With Backoff"); xlabel("Input Amplitude"); ylabel("Phase Modulation (Deg)");
  grid on; set(gca,"fontsize",20);
endif

% Calculate EVM
x_evm = x_bo(1+wola_len/2:end);
y_evm = y_bo(1+wola_len/2:end);
cfg_evm.en_plot = en_plot;
evm_bo = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
snr_bo = -20*log10(evm_bo/100);

% Plot PSDs
if en_plot
  px_bo = calculate_psd(x_bo,fs,rbw); px_bo = scale_psd(px_bo,f,bw,scs,start_sc,num_sc);
  py_bo = calculate_psd(y_bo,fs,rbw); py_bo = scale_psd(py_bo,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px_bo),'linewidth',5); plot(f,10*log10(py),'linewidth',3); plot(f,10*log10(py_bo),'linewidth',1);
  title("PA PSD With Backoff"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

% Calculate PAPR
[papr_x_bo,avg_x_bo,peak_x_bo] = calculate_papr(x_bo,100);
[papr_y_bo,avg_y_bo,peak_y_bo] = calculate_papr(y_bo,100);

%{
disp("Input power without power backoff minus input power with power backoff");
disp(10*log10(avg_x/avg_x_bo));
disp("Input peak power without power backoff minus input peak power with power backoff");
disp(10*log10(peak_x/peak_x_bo));
disp("PA power without power backoff minus PA power with power backoff");
pin_backoff = avg_y/avg_y_bo; % linear power backoff (for next experiment)
disp(10*log10(pin_backoff));
disp("PA peak power without power backoff minus PA peak power with power backoff");
disp(10*log10(peak_y/peak_y_bo));
%}


% Output metrics
disp("--- PA EVM / With DPD / With Backoff ---");
disp(evm); disp(evm_dpd); disp(evm_bo); disp("");

disp("--- PA SNR / With DPD / With Backoff ---");
disp(snr); disp(snr_dpd); disp(snr_bo); disp("");

disp("--- PAPR: Input / Output ---");
disp(papr_x); disp(papr_y); disp("");

disp("--- PAPR With DPD: Input / Output ---");
disp(papr_x_dpd); disp(papr_y_dpd); disp("");

disp("--- PAPR With Backoff: Input / Output ---");
disp(papr_x_bo); disp(papr_y_bo); disp("");

disp("--- Average PA Power / With DPD / With Backoff ---");
disp(10*log10(avg_y)); disp(10*log10(avg_y_dpd)); disp(10*log10(avg_y_bo)); disp("");

disp("--- Peak PA Power / With DPD / With Backoff ---");
disp(10*log10(peak_y)); disp(10*log10(peak_y_dpd)); disp(10*log10(peak_y_bo)); disp("");
