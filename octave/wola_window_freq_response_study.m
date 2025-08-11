% Plot frequency response of WOLA window

clear; clc; close all;

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




sym_len = 2048+144;
sym_len = 2048;
wola_lens = 0:16:144;
nfft = sym_len*2^5;

figure; hold on;
for wdx = 1:length(wola_lens)
  wola_len = wola_lens(wdx);
  b = generate_wola_window(wola_len,sym_len);
  b = b/sum(b);
  [h,w] = freqz(b,1,nfft);
  plot(w,20*log10(abs(h)),'linewidth',5);
endfor
title("WOLA Window Magnitude Response");
xlabel("Digital Frequency"); ylabel("Magnitude Response");
lobj = legend(mat2str(wola_lens));
set(lobj,'fontsize',25);
xlim([0 0.25]); ylim([-75 0]); grid on;
set(gca,"fontsize",25);