% Frequency domain representation of DT impulse train

clear; clc; close all;

N = 5; % periodicity

x = zeros(1,N); x(1) = 1;
x = [x x];

x_fft = fft(x);
%f = 0:2/N:(2-2/N);
f = 0:2/length(x):(2-2/length(x));

figure; stem(f,abs(x_fft));