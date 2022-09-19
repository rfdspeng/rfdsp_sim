% Author: Ryan Tsai
% Created: 2022-05-06

function [y,cfg] = resample_custom(x,p,q,nbw)
% Custom resample function

% Unless otherwise specified, assume normalized signal bandwidth prior to resampling is <= 95%
if ~exist('nbw','var'), nbw = 0.95; endif

% Estimate number of taps for interpolation/AA filter
r = max(p,q); % determines the stopband
transition = [nbw 2-nbw]/r;
ripple = 0.001; atten = 120; % dB specs
dev = [1-10^(-ripple/20) 10^(-atten/20)];
n = kaiserord(transition,[1 0],dev);
n_fh = ceil(atten/22/diff(transition)); % fred harris

% Generate filter
n = n+mod(n,2);
weight = max(dev)./dev;
weight = weight.*weight;
b = firls(n,[0 transition 1],[1 1 0 0],weight); b = b/sum(b);

% Filter sanity check
[h,w] = freqz(b,1,2^12); h = abs(h); w = w/pi;
hpb = h(w <= transition(1)); hsb = h(w >= transition(2));
max_ripple = max(abs(20*log10(hpb)));
min_atten = min(-20*log10(hsb));
if max_ripple > ripple || min_atten < atten, error("resample_custom filter does not meet spec"); endif

% Resample
y = upsample(x,p);
y = filter(b,1,[y zeros(1,n/2)]);
y = y(1+n/2:end);
y = downsample(y,q);

% Debug
cfg.n = n;
cfg.f = [0 transition 1];
cfg.a = [1 1 0 0];
cfg.w = weight;
cfg.max_ripple = max_ripple;
cfg.min_atten = min_atten;

endfunction
