% Author: Ryan Tsai
% Created: 2022-03-30

function [x,x_standard,cfg] = ofdm_wavgen(nsym,bw,scs,num_sc,start_sc,modorder,en_tprecode,ncp,wola,seed)
% Generates OFDM waveform with 1-subcarrier resolution in subcarrier allocation
% bw = BW in MHz
% scs = SCS in kHz
% ncp = CP length in percentage of nfft
% wola = WOLA length in percentage of nfft

if ~exist('ncp','var'), ncp = 0; endif
if ~exist('wola','var'), wola = 0; endif
if exist('seed','var'), rand("state",seed); endif

nrb = bw*5; % max number of resource blocks
nfft = 2^ceil(log2(nrb*12));
ncp = round(nfft*ncp/100);
sym_len = nfft+ncp;
fs = nfft*scs/1000; % waveform sampling rate
bps = log2(modorder); % bits per tone
wola_len = round(nfft*wola/100);
wola_len = wola_len+mod(wola_len,2);
b = generate_wola_window(wola_len,sym_len);

x_standard = zeros(1,nsym*sym_len); % no WOLA
x = zeros(1,nsym*sym_len+wola_len); % with WOLA
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
  x_standard(1+(sdx-1)*sym_len:sdx*sym_len) = sym_cp;
  
  sym_wola = [sym(nfft-ncp-wola_len/2+1:end) sym sym(1:wola_len/2)];
  sym_wola = sym_wola.*b;
  x(1+(sdx-1)*sym_len:sdx*sym_len+wola_len) = sym_wola + x(1+(sdx-1)*sym_len:sdx*sym_len+wola_len);
  
  if sdx == 1
    tone_idx = ones(1,nrb*12);
    tone_idx(1+start_sc:start_sc+num_sc) = 2;
    tone_idx = [zeros(1,(nfft-nrb*12)/2) tone_idx zeros(1,(nfft-nrb*12)/2)];
    tone_idx = [tone_idx(nfft/2+1:end) tone_idx(1:nfft/2)];
  endif
endfor

%x1 = [zeros(1,wola_len/2) x_standard zeros(1,wola_len/2)];
%len = 2*(wola_len+nfft);
%figure; hold on; plot(abs(x1(1:len))); plot(abs(x(1:len)));

% Struct for calculating EVM
cfg.fs = fs; cfg.nfft = nfft; cfg.ncp = ncp; cfg.nrb = nrb; cfg.nsym = nsym;
cfg.bw = bw; cfg.scs = scs; cfg.num_sc = num_sc; cfg.start_sc = start_sc;
cfg.modorder = modorder; cfg.en_tprecode = en_tprecode; cfg.tone_idx = tone_idx;
cfg.wola_len = wola_len;

endfunction

function tones = modulation_mapper(bitstream,modorder)
  if modorder == 4
    tones = 1-2*bitstream(1:2:end) + 1j*(1-2*bitstream(2:2:end));
    tones = tones/sqrt(2);
  endif
endfunction

function b = generate_wola_window(wola_len,sym_len)
% sym_len = nfft+ncp

if wola_len > 0
  b = hanning(2*wola_len).';
  b = [b(1:wola_len) ones(1,sym_len-wola_len) b(wola_len+1:end)];
else
  b = ones(1,sym_len);
endif
%b = hanning(wola_len+sym_len).';
  
endfunction
