% Author: Ryan Tsai
% Created: 2022-02-06

function [p,f] = calculate_psd(x,fs,rbw)
nfft = ceil(fs/rbw);
[p,f] = pwelch(x,hanning(nfft),[],nfft,fs,'whole');
p = fftshift(p); f = f-fs/2;
endfunction
