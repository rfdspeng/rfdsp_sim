% Evaluate WOLA impact on EVM
% Evaluate window impact on a symbol
% Evaluate ISI impact on a symbol
% Evaluate both effects

clear; clc; close all;

nsym = 50; bw = 20; scs = 15; num_sc = 100; start_sc = 400; modorder = 4; en_tprecode = 0;
ncp = 0;
wolas = (0:16:144)/2048*100;
evms = zeros(size(wolas));
for wdx = 1:length(wolas)
  wola = wolas(wdx);
  [x,x_standard,cfg_evm] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola);
  x_evm = x_standard;
  y_evm = x(cfg_evm.wola_len/2+1:end);
  
  if wdx == 1 || wdx == length(wolas), cfg_evm.en_plot = 1;
  else, cfg_evm.en_plot = 0; endif
  evm = ofdm_evm_calculator(cfg_evm,x_evm,y_evm);
  evms(wdx) = evm;
endfor