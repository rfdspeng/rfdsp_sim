% Author: Ryan Tsai
% Created: 2022-04-02

function evm = ofdm_evm_calculator(cfg,x,y)
% Calculates EVM (%) for OFDM waveform
% x = reference, y = signal
% Signals must be at baseband sampling rate and time-aligned

if ~isfield(cfg,'en_plot'), en_plot = 0;
else, en_plot = cfg.en_plot; endif

n = min(length(x),length(y));
nfft = cfg.nfft; ncp = cfg.ncp; sym_len = nfft+ncp;
en_tprecode = cfg.en_tprecode;
fft_start = max(round(ncp/2),1);
tones_cell = cell(2,1); symdx = 1;
tone_idx = cfg.tone_idx;
tone_idx = tone_idx == 2;
total_signal_power = 0;
total_error_power = 0;
for sdx = 1:sym_len:n-sym_len+1
  sym_x = x(sdx:sdx+sym_len-1);
  sym_y = y(sdx:sdx+sym_len-1);
  
  yy = sym_y.'; xx = sym_x.';
  xx_hat = (xx'*yy)/(yy'*yy)*yy;
  sym_y = xx_hat.';
  
  sym_x = sym_x(fft_start:fft_start+nfft-1);
  sym_y = sym_y(fft_start:fft_start+nfft-1);
  nshift = ncp-fft_start+1;
  sym_x = circshift(sym_x,-nshift);
  sym_y = circshift(sym_y,-nshift);
  
  tones_x = fft(sym_x); tones_y = fft(sym_y);
  tones_x = tones_x(tone_idx);
  tones_y = tones_y(tone_idx); 
  
  if en_tprecode
    tones_x = ifft(tones_x);
    tones_y = ifft(tones_y);
  endif
  
  signal_power = tones_x*tones_x';
  err = tones_x - tones_y;
  error_power = err*err';
  total_signal_power = total_signal_power+signal_power;
  total_error_power = total_error_power+error_power;
  
  tones_cell{1,symdx} = tones_x;
  tones_cell{2,symdx} = tones_y;
  symdx = symdx+1;
endfor
evm = sqrt(total_error_power/total_signal_power)*100;

if en_plot
  tones_x = [tones_cell{1,:}];
  tones_y = [tones_cell{2,:}];
  a_norm = rms(tones_x);
  tones_x = tones_x/a_norm;
  tones_y = tones_y/a_norm;
  figure; hold on;
  plot(tones_x,'x','markersize',25);
  plot(tones_y,'x','markersize',10);
  title("Constellation");
  xlabel("I"); ylabel("Q");
  set(gca,"fontsize",25);
endif

endfunction
