% How to simulate two-tone test using complex baseband equivalents?

clear; clc; close all;

en_plot = 1;
niter = 1;

xlen = 2^17; nfft = round(xlen/2^6);
% The results are wonky if frf = pi/2, I think because the 3rd harmonic term wraps around and lands on IM3
%fbb = 2*pi/40;
%frf = pi/2;
fbb = 2*pi/160;
fbb = 2*pi/4;
frf = pi/8;

piip3 = 0;
pin = -40;
a1 = 1;

r = 10^(pin/20);
a3 = -4/3*a1*10^(-piip3/10);



if niter > 1, en_plot = 0; endif
noises = zeros(1,niter);
noises_rf = zeros(1,niter);
for iterdx = 1:niter
  disp(iterdx);
  xbb = exp(1j*(fbb*(1:xlen)+rand()*2*pi)) + exp(1j*(-fbb*(1:xlen)+rand()*2*pi));

  hannwin = hanning(length(xbb)); hannwin = hannwin.';
  xbbwin = xbb.*hannwin;
  xbb_psd = 20*log10(abs(fftshift(fft(xbbwin)))); xbb_psd = xbb_psd-max(xbb_psd);
  fbb_psd = linspace(-1,1,length(xbb_psd));
  if en_plot, figure; plot(fbb_psd,xbb_psd); grid on; endif

  [xbb_psd2,fbb_psd2] = pwelch(xbb,hanning(nfft),[],nfft,2*pi);
  xbb_psd2 = fftshift(xbb_psd2); fbb_psd2 = fbb_psd2-pi;
  xbb_psd2 = 10*log10(xbb_psd2); xbb_psd2 = xbb_psd2-max(xbb_psd2);
  figure; plot(fbb_psd2,xbb_psd2); grid on;


  lo = exp(1j*(frf*(1:xlen)+rand()*2*pi));
  xrf = real(xbb.*lo);
  xrfwin = xrf.*hannwin;
  xrf_psd = 20*log10(abs(fftshift(fft(xrfwin)))); xrf_psd = xrf_psd-max(xrf_psd);
  frf_psd = linspace(-1,1,length(xrf_psd));
  if en_plot, figure; plot(frf_psd,xrf_psd); grid on; endif

  %[xrf_psd1,frf_psd1] = pwelch_custom(xrf,hanning(nfft));
  [xrf_psd1,frf_psd1] = pwelch_custom(xrf,ones(1,nfft));
  xrf_psd1 = 10*log10(xrf_psd1); xrf_psd1 = xrf_psd1-max(xrf_psd1);
  figure; plot(frf_psd1,xrf_psd1); grid on;

  [xrf_psd2,frf_psd2] = pwelch(xrf,hanning(nfft),[],nfft,2*pi);
  xrf_psd2 = fftshift(xrf_psd2); frf_psd2 = frf_psd2-pi;
  xrf_psd2 = 10*log10(xrf_psd2); xrf_psd2 = xrf_psd2-max(xrf_psd2);
  figure; plot(frf_psd2,xrf_psd2); grid on;




  yrf = a1*r*xrf + a3*r^3*xrf.^3;
  yrfwin = yrf.*hannwin;
  yrf_psd_lin = abs(fftshift(fft(yrfwin))).^2;
  yrf_psd = 10*log10(yrf_psd_lin); yrf_psd = yrf_psd-max(yrf_psd);
  if en_plot, figure; plot(frf_psd,yrf_psd); grid on; endif


  ybb = a1*r*xbb + 3/4*a3*r^3*abs(xbb).^2.*xbb;
  ybbwin = ybb.*hannwin;
  ybb_psd_lin = abs(fftshift(fft(ybbwin))).^2;
  ybb_psd = 10*log10(ybb_psd_lin); ybb_psd = ybb_psd-max(ybb_psd);
  if en_plot, figure; plot(fbb_psd,ybb_psd); grid on; endif


  sl = -fbb-fbb/4; sh = -fbb+fbb/4;
  nl = -3*fbb-fbb/4; nh = -3*fbb+fbb/4;
  sl = sl/pi; sh = sh/pi; nl = nl/pi; nh = nh/pi;
  noise_dbc = calculate_noise_dbc(ybb_psd_lin,fbb_psd,sl,sh,nl,nh,en_plot);
  noises(iterdx) = noise_dbc;


  sl = +frf-fbb-fbb/4; sh = +frf-fbb+fbb/4;
  nl = +frf-3*fbb-fbb/4; nh = +frf-3*fbb+fbb/4;
  sl = sl/pi; sh = sh/pi; nl = nl/pi; nh = nh/pi;
  noise_dbc_rf = calculate_noise_dbc(yrf_psd_lin,frf_psd,sl,sh,nl,nh,en_plot);
  noises_rf(iterdx) = noise_dbc_rf;
endfor

noise_avg = mean(noises)
noise_avg_rf = mean(noises_rf)
