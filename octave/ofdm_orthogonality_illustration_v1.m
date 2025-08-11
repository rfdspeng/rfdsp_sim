% Illustrate how subcarriers are orthogonal in OFDM
% Of course, this is only true when you apply a rectangular window of length nfft to each OFDM symbol,
% where nfft is the FFT size used to generate the OFDM symbols

clear; clc; close all;

nsc = 5; % number of subcarriers to simulate

nsym = 140*10; bw = 20; scs = 15; num_sc = 1; start_sc = 41*12; modorder = 4; en_tprecode = 0;
ncp = 0; wola = 5;

psds = cell(1,nsc);
psds_wola = cell(1,nsc);
for sdx = 0:nsc-1
  [x,x_standard,cfg] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc+sdx,modorder,en_tprecode,ncp,wola);
  nfft = cfg.nfft;
  fs = 2;
  [p,f] = pwelch_custom(x_standard,hanning(nfft*16),[],nfft*16,fs);
  p = p/max(p);
  psds{sdx+1} = p;
  [p,f] = pwelch_custom(x,hanning(nfft*16),[],nfft*16,fs);
  p = p/max(p);
  psds_wola{sdx+1} = p;
endfor

figure; hold on;
for pdx = 1:length(psds)
  plot(f,psds{pdx});
endfor
%xlim([-fs/2 fs/2]); grid on;
xlim([-0.12 -0.09]);

figure; hold on;
for pdx = 1:length(psds_wola)
  plot(f,psds_wola{pdx});
endfor
%xlim([-fs/2 fs/2]); grid on;
xlim([-0.12 -0.09]);

figure; hold on;
for pdx = 1:length(psds)
  plot(f,10*log10(psds{pdx}));
endfor
%xlim([-fs/2 fs/2]); grid on;
xlim([-0.5 0.5]);
ylim([-100 0]);
grid on;

figure; hold on;
for pdx = 1:length(psds_wola)
  plot(f,10*log10(psds_wola{pdx}));
endfor
%xlim([-fs/2 fs/2]); grid on;
xlim([-0.5 0.5]);
ylim([-100 0]);
grid on;