% Sanity check ofdm_wavgen() and ofdm_evm_calculator()

clear; clc; close all;

en_plot = 1;

nsym = 30; bw = 20; scs = 15; num_sc = 100; start_sc = 400; modorder = 4; en_tprecode = 0;
ncp = 0; wola = 5; wola = 0;

[x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola);
x_evm = x_standard;
y_evm = x(cfg_evm.wola_len/2+1:end);

cfg_evm.en_plot = en_plot;
evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
disp(evm);