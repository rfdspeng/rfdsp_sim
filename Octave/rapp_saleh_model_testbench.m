clear; clc; close all;

% Rapp params
cfg.g = 10^(30/20);
cfg.smoothness = 2;
cfg.osat = 25;

% Saleh params
cfg.a = 0.5; cfg.b = 10;

x = linspace(0,1,2^12);
y = rapp_saleh_model(cfg,x);
figure; plot(abs(x),abs(y));
figure; plot(abs(x),20*log10(abs(y)./abs(x)));
figure; plot(abs(x),unwrap(angle(y))*180/pi);

len = 2^12; fs = 30.72; bw = 5;
cfg_wav.fs = fs; cfg_wav.bw = bw; cfg_wav.length = len;
x = wavgen(cfg_wav); x = x/max(abs(x));
y = rapp_saleh_model(cfg,x);

nfft = round(length(x)/8);
[p,f] = pwelch(x,hanning(nfft),[],nfft,fs);
py = pwelch(y,hanning(nfft),[],nfft,fs);
p = fftshift(p); py = fftshift(py); f = f-fs/2;
p = p/max(p); py = py/max(py);
figure; hold on;
plot(f,10*log10(p),'LineWidth',10);
plot(f,10*log10(py),'LineWidth',5);
xlim([-fs/2 fs/2]); ylim([-100 0]);
yticks(-75:5:0);
%xline(-bw/2); xline(bw/2);
%xline(-bw-bw/2); xline(bw+bw/2);
grid on;

mbw = 4.5;
psig = py(f >= -bw/2 & f <= bw/2);
pnoise_low = py(f >= -bw-mbw/2 & f <= -bw+mbw/2);
pnoise_high = py(f >= bw-mbw/2 & f <= bw+mbw/2);
psig = 10*log10(sum(psig)); pnoise_low = 10*log10(sum(pnoise_low)); pnoise_high = 10*log10(sum(pnoise_high));
nlowdb = pnoise_low-psig
nhighdb = pnoise_high-psig