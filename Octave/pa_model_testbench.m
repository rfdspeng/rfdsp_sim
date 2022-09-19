% Test bench for generating a "realistic" PA model

% Assume APT mode needs to support up to 17dBm without saturating for a max PA bias of 3.5V

clear; clc; close all;

addpath("tools");
addpath("models");

en_plot = 1;

% Waveform params
en_tprecode = 1;
target_papr = 4.5;
osr = 4; % oversampling ratio after waveform generation

% PA params
gain = 30; % dB
pmax_apt = 17; % RF, dBm


v = sqrt(2)*10^((pmax_apt-gain)/20); % complex baseband rms voltage @ PA input

% Waveform
nsym = 120; bw = 20; scs = 15; num_sc = 600; start_sc = 600-num_sc/2; modorder = 4;
ncp = 10; wola = 2; seed = 1;
[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
wola_len = cfg_evm.wola_len;

% Upsample waveform prior to CFR
x = resample_custom(x,osr,1);
fs = fs*osr;

% Clip waveform
x_cfr = cfr(x,target_papr);

% Scale waveform to PA input
x_cfr = x_cfr/rms(x_cfr)*v;


% PA model
cfg_pa.g = 10^(gain/20);
cfg_pa.smoothness = 2;
cfg_pa.osat = 36;
cfg_pa.a = 0.5; cfg_pa.b = 10;
y = rapp_saleh_model(cfg_pa,x_cfr);

papr = calculate_papr(x,99.99);
papr_clip = calculate_papr(x_cfr,99.99);
papr_pa = calculate_papr(y,99.99);
comp = calculate_compression(x_cfr,y,en_plot);

disp("--- PAPR (dB): Unclipped / Clipped / PA Output ---");
disp(papr); disp(papr_clip); disp(papr_pa); disp("");

disp("--- Peak Compression @ PA Output (dB) ---");
disp(comp); disp("");

% Plot PSDs
if en_plot
  rbw = scs/1000/2^6;
  [px,f] = calculate_psd(x,fs,rbw); px = scale_psd(px,f,bw,scs,start_sc,num_sc);
  pcfr = calculate_psd(x_cfr,fs,rbw); pcfr = scale_psd(pcfr,f,bw,scs,start_sc,num_sc);
  py = calculate_psd(y,fs,rbw); py = scale_psd(py,f,bw,scs,start_sc,num_sc);
  figure; hold on;
  plot(f,10*log10(px),'linewidth',5);
  plot(f,10*log10(pcfr),'linewidth',3);
  plot(f,10*log10(py),'linewidth',1);
  title("PSD"); xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
  legend("Mod","CFR","PA",'location','south');
  xlim([-fs/2 fs/2]);
  grid on; set(gca,"fontsize",20);
endif

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
