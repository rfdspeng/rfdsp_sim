% Confirm that SQNR simulation matches the equations

clear; clc; close all;

en_plot = 1;

bitwidth = 16; setpoint = -16; % data params

nsym = 120; bw = 20; scs = 15; num_sc = 100; start_sc = 600-num_sc/2; modorder = 4; en_tprecode = 0;
ncp = 10; wola = 2; %wola = 0;

[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola);
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
disp(snr);
%disp(evm_td);
disp(snr_td);

snr_expected = 6.02*bitwidth+1.76+(setpoint+0)+10*log10(fs/(num_sc*scs/1000));
snr_td_expected = 6.02*bitwidth+1.76+(setpoint+0);
disp(snr_expected); disp(snr_td_expected);

rbw = scs/1000/2^3;
[p,f] = calculate_psd(x,fs,rbw);
p_fl = calculate_psd(x_fl,fs,rbw);
p = scale_psd (p,f,bw,scs,start_sc,num_sc);
p_fl = scale_psd (p_fl,f,bw,scs,start_sc,num_sc);
figure; hold on; plot(f,10*log10(p)); plot(f,10*log10(p_fl));
set(gca,"ytick",-120:2.5:5,"fontsize",20);
xlim([-fs/2 fs/2]); ylim([-120 5]);
grid on;