clear; clc; close all;

% Generate input
nbw = 0.05;
len = 2^17; fs = 2; rbw = 2/2^12;
x1 = randn(1,len) + 1j*randn(1,len);
b = firls(300,[0 nbw nbw+0.05 1],[1 1 0 0]);
figure; freqz(b);
x = filter(b,1,x1);
[p,f] = calculate_psd(x,fs,rbw);
figure; plot(f,10*log10(p));

% Generate blocker
nbw_blocker = 0.2;
x1 = randn(1,len) + 1j*randn(1,len);
b = firls(300,[0 nbw_blocker nbw_blocker+0.05 1],[1 1 0 0]);
figure; freqz(b);
x_blocker = filter(b,1,x1);
x_blocker = x_blocker/rms(x_blocker)*rms(x)*10;
x_blocker = x_blocker.*exp(1j*2*pi/4*(1:len));
[p,f] = calculate_psd(x_blocker,fs,rbw);
figure; plot(f,10*log10(p));

% Add signal+blocker
x = x+x_blocker;
[p,f] = calculate_psd(x,fs,rbw);
figure; plot(f,10*log10(p));

% Antialiasing filter design
wp = nbw; ws = 2/4-nbw_blocker;
ripple_spec = 0.1; % dB
rej_spec = 50; % dB
dev = [1-10^(-ripple_spec/20) 10^(-rej_spec/20)];
w = max(dev)./dev;
w = w.*w;
b = firls(50,[0 wp ws 1],[1 1 0 0],w);
figure; freqz(b);
y = filter(b,1,x);
y = downsample(y,4);
[p,f] = calculate_psd(y,fs,rbw);
figure; plot(f,10*log10(p));


function [p,f] = calculate_psd(x,fs,rbw)
nfft = ceil(fs/rbw);
[p,f] = pwelch(x,hanning(nfft),[],nfft,fs);
p = fftshift(p); f = f-fs/2;
endfunction