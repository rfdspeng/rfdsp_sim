% Author: Ryan Tsai
% Created: 2022-03-30

function [x,fs,nfft,ncp,nrb] = ofdm_wavgen_by_subcarrier(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,seed)
% Generates OFDM waveform with 1-subcarrier resolution in subcarrier allocation
% bw = BW in MHz
% scs = SCS in kHz

nrb = bw*5; % max number of resource blocks
nfft = 2^ceil(log2(nrb*12));
ncp = nfft*288/4096; ncp = 0; %ncp = nfft;
sym_len = nfft+ncp;
fs = nfft*scs/1000; % waveform sampling rate
bps = log2(modorder); % bits per tone

if exist('seed','var'), rand("state",seed); endif

x = zeros(1,nsym*sym_len);
for sdx = 1:nsym
  bitstream = rand(1,bps*num_sc)-0.5>0;
  tones = modulation_mapper(bitstream,modorder);
  if en_tprecode, tones = fft(tones); endif
  tones_nrb = zeros(1,nrb*12);
  tones_nrb(1+start_sc:start_sc+num_sc) = tones;
  tones_all = [zeros(1,(nfft-nrb*12)/2) tones_nrb zeros(1,(nfft-nrb*12)/2)];
  tones_all = [tones_all(nfft/2+1:end) tones_all(1:nfft/2)];
  sym = ifft(tones_all);
  sym_cp = [sym(nfft-ncp+1:end) sym];
  x(1+(sdx-1)*sym_len:sdx*sym_len) = sym_cp;
endfor

endfunction

function tones = modulation_mapper(bitstream,modorder)
  if modorder == 4
    tones = 1-2*bitstream(1:2:end) + 1j*(1-2*bitstream(2:2:end));
    tones = tones/sqrt(2);
  endif
endfunction
