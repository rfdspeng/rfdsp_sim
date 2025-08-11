% Author: Ryan Tsai
% Created: 2022-03-30

function [p,f,num_segments] = pwelch_custom(x,win,overlap_perc,nfft,fs,center)
% Octave built-in pwelch shows high DC offset for some reason
% This is a custom function that implements pwelch
% x = signal
% win = window. Window length <= FFT size
% overlap_perc = overlap percentage in %. Default 50%.
% nfft = FFT size. Default = window length.
% fs = sampling rate. Default = 2*pi.
% center: If 'centered', return PSD centered at DC

if ~exist('nfft','var'), nfft = length(win); endif
if length(win) > nfft, error("Window length > FFT size"); endif
if isempty(overlap_perc), overlap_perc = 50; endif
if ~exist('fs','var'), fs = 2*pi; endif
if ~exist('center','var'), center = 'centered'; endif
if iscolumn(x), x = x.'; endif
if iscolumn(win), win = win.'; endif

nsegment = length(win);
overlap_perc = overlap_perc/100;
step_size = round(nsegment*overlap_perc);
p = zeros(1,nfft); % PSD in V^2/bin
f = 0:fs/nfft:(fs-fs/nfft);

num_segments = 0;
for sdx = 1:step_size:(length(x)-nsegment)
  xwin = x(sdx:sdx+nsegment-1).*win;
  xwin = [xwin zeros(1,nfft-nsegment)];
  xfft = fft(xwin);
  xpsd = xfft.*conj(xfft);
  p = p+xpsd;
  num_segments = num_segments+1;
endfor
p = p/num_segments;

if center == 'centered'
  p = fftshift(p);
  f = f-fs/2;
endif

endfunction
