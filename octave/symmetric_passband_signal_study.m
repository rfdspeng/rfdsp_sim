% Symmetric passband signal (amplitude modulation only) upconverted to RF
% Sample directly at RF frequency

clear; clc; close all;

nsym = 2^9; bps = 3; osr = 21; osr = 2;
x = mod_ask(nsym,bps,osr);
figure; plot(x(1:100));
hannwin = hanning(length(x)).';
xwin = hannwin.*x;
xwin = x;
x_psd = abs(fftshift(fft(xwin))).^2;
figure; plot(10*log10(x_psd));