% Author: Ryan Tsai
% Created: 2022-01-29

function [y,outstruct] = tow_thomas_lpf(cfg,x)
% https://en.wikipedia.org/wiki/Electronic_filter_topology%Biquad_filter_topology
% cfg.g_dc: DC gain (linear)
% cfg.f0: Natural frequency
% cfg.Q: Quality factor
% cfg.fs: Sampling rate of x
% cfg.rbw: RBW for returning frequency response
% cfg.opt: "impinvar" or "bilinear"

if ~isfield(cfg,'opt'), cfg.opt = "impinvar"; endif
g_dc = cfg.g_dc; f0 = cfg.f0; Q = cfg.Q; fs = cfg.fs; rbw = cfg.rbw; opt = cfg.opt;
nfft = ceil(fs/2/rbw); % 0 to fs/2

% Analog coefficients
w0 = 2*pi*f0;
b = g_dc*w0^2;
a = [1 w0/Q w0^2];

% CT-DT transform
if strcmp(opt,"impinvar"), [bz,az] = impinvar(b,a,fs);
elseif strcmp(opt,"bilinear"), [bz,az] = bilinear(b,a); % probably not correct
endif
bz = bz/sum(bz)*sum(az)*g_dc;

% Analog frequency response
w = 2*pi*rbw*(0:nfft-1); f = w/2/pi;
h = freqs(b,a,w);
h_theory = g_dc*w0^2./( (1j*w).^2 + w0/Q*1j*w + w0^2 );

% Digital frequency response
[hz,wz] = freqz(bz,az,nfft);

figure; hold on;
plot(f,20*log10(abs(h)));
plot(f,20*log10(abs(h_theory)));
plot(f,20*log10(abs(hz)));

outstruct.b = b; outstruct.a = a;
outstruct.bz = bz; outstruct.az = az;

y = filter(bz,az,x);

endfunction
