% Plot frequency response of WOLA window

clear; clc; close all;

nsym = 14*40; bw = 20; scs = 15; num_sc = 800; start_sc = 200; modorder = 4; en_tprecode = 0;
ncp = 144/2048*100; ncp = 0;

rbw = scs/1000/2^6;
wolas = (0:16:144)/2048*100;
figure; hold on;
for wdx = 1:length(wolas)
  wola = wolas(wdx);
  [x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola);
  fs = cfg_evm.fs;
  disp(cfg_evm.wola_len);
  [p,f] = calculate_psd(x,fs,rbw);
  p = scale_psd (p,f,bw,scs,start_sc,num_sc);
  plot(f,10*log10(p));
endfor
title("Spectral Leakage vs. WOLA Length");
xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
%legend(mat2str(wola_lens));
xlim([-fs/2 fs/2]); ylim([-125 25]); grid on;
set(gca,"fontsize",25);