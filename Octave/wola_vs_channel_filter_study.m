% Compare WOLA and channel filter performance

clear; clc; close all;

en_plot = 1;

nsym = 14*5; bw = 20; scs = 15; num_sc = 800; start_sc = 200; modorder = 4; en_tprecode = 0;
ncp = 0; wola = 64/2048*100; seed = 1;

rbw = scs/1000/2^6;

% WOLA
[x_wola,x,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed);
fs = cfg_evm.fs;
x_evm = x;
y_evm = x_wola(cfg_evm.wola_len/2+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
disp(evm);

% Compare standard to WOLA PSD
[p,f] = calculate_psd(x,fs,rbw); p = scale_psd(p,f,bw,scs,start_sc,num_sc);
p_wola = calculate_psd(x_wola,fs,rbw); p_wola = scale_psd(p_wola,f,bw,scs,start_sc,num_sc);
figure; hold on;
plot(f,10*log10(p)); plot(f,10*log10(p_wola));
title("Spectral Leakage vs. WOLA Length");
xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
%legend(mat2str(wola_lens));
xlim([-fs/2 fs/2]); ylim([-125 10]); grid on;
%set(gca,"fontsize",25);
set(gca,"ytick",-125:5:10);
set(gca,"xtick",-15:1:15);

% Define channel filter requirements
fpass = num_sc*scs/1000/fs; % passband
% Stopband requirements: @f_stop_i, attenuation must be a_db_i (in dB)
fstop = (14:2:30)/30.72;
a_db = [-35 -47 -53 -59 -62 -65 -68 -70 -73];

% Generate channel filter coefficients
%dev = [1-10^(-0.01/20) 10.^(a_db/20)];
%weight = max(dev)./dev; weight = weight.*weight;
a = 10.^(a_db/20);
fstop_interp = linspace(min(fstop),max(fstop),2^12-1);
a_interp = interp1(fstop,a,fstop_interp);
n = 120; f = [0 fpass fstop_interp 1]; a = [1 1 a_interp 0];
b = firls(n,f,a); %b = b/sum(b);
[h,w] = freqz(b,1,2^11*2^5); w = w/pi;
hpass = h(w <= fpass);
ripple = max(abs(20*log10(abs(hpass))));
disp(ripple);
figure; plot(w,20*log10(abs(h))); grid on;
set(gca,"xtick",0:0.05:1);

% Run standard waveform through channel filter
gd = n/2;
x1 = [x zeros(1,length(b))];
x_ch = filter(b,1,x1);
[p_ch,f] = calculate_psd(x_ch,fs,rbw); p_ch = scale_psd(p_ch,f,bw,scs,start_sc,num_sc);
figure; hold on;
plot(f,10*log10(p)); plot(f,10*log10(p_ch)); plot(f,10*log10(p_wola));
title("Spectral Leakage");
xlabel("Frequency (MHz)"); ylabel("PSD (dB)");
legend("Default Waveform","Channel Filter","WOLA");
xlim([-fs/2 fs/2]); ylim([-125 10]); grid on;
%set(lobj,"fontsize",10);
set(gca,"fontsize",25);
%set(gca,"ytick",-125:5:10);
%set(gca,"xtick",-15:1:15);

% Plot time domain samples
%figure; hold on;
%plot([zeros(1,gd) abs(x)]);
%plot(abs(x_ch));

% Calculate EVM for channel filter
x_evm = x;
y_evm = x_ch(gd+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
disp(evm);