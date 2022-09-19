% Illustrate how subcarriers are orthogonal in OFDM
% Of course, this is only true when you apply a rectangular window of length nfft to each OFDM symbol,
% where nfft is the FFT size used to generate the OFDM symbols

clear; clc; close all;

nfft = 4;
nfft_plot = nfft*2^6;
w = 0:2*pi/nfft_plot:2*pi; w = w(1:end-1);

subcarriers = zeros(nfft,nfft);

figure; hold on;
for ndx = 1:nfft
  subcarrier = exp(1j*(ndx-1)*2*pi/nfft*(0:nfft-1));
  subcarriers(ndx,:) = subcarrier;
  sfft = fft([subcarrier zeros(1,nfft_plot-nfft)]);
  plot(w,abs(sfft)/nfft,'linewidth',7);
endfor
title("Subcarriers in Frequency Domain");
xlabel("Digital Frequency"); ylabel("|FFT|");
xlim([min(w) max(w)]); grid on;
set(gca,"fontsize",25);
%set(gca,'linewidth',5);

% Rectangular window
figure;
subcarrier = ones(1,nfft);
sfft = fft([subcarrier zeros(1,nfft_plot-nfft)]);
plot(w,abs(sfft)/nfft,'linewidth',7);
title("Rectangular Window");
xlabel("Digital Frequency"); ylabel("Magnitude Response");
xlim([min(w) max(w)]); grid on;
set(gca,"fontsize",25);
%set(gca,'linewidth',5);