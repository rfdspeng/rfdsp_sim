% Author: Ryan Tsai
% Created: 2022-05-06

function [b,cfg] = firls_wrapper(n,fpass,fstop,ripple,atten,en_check)
% Wrapper function to auto-generate a linear-phase LPF
% fpass = normalized passband frequency
% fstop = normalized stopband frequency
% ripple = ripple spec in dB
% atten = attenuation spec in dB
% en_check = sanity check that specs are met

if ~exist('en_check','var'), en_check = 1; endif

dev = [1-10^(-ripple/20) 10^(-atten/20)];
weight = max(dev)./dev;
weight = weight.*weight;
if isempty(n)
  n = kaiserord([fpass fstop],[1 0],dev);
  n = n+mod(n,2); % even-order filter
endif
b = firls(n,[0 fpass fstop 1],[1 1 0 0],weight); b = b/sum(b);

% Debug
cfg.n = n;
cfg.f = [0 fpass fstop 1];
cfg.a = [1 1 0 0];
cfg.w = weight;

if en_check
  [h,w] = freqz(b,1,2^12); h = abs(h); w = w/pi;
  hpb = h(w <= fpass); hsb = h(w >= fstop);
  max_ripple = max(abs(20*log10(hpb)));
  min_atten = min(-20*log10(hsb));
  if max_ripple > ripple || min_atten < atten, error("firls_wrapper filter does not meet spec"); endif

  cfg.max_ripple = max_ripple;
  cfg.min_atten = min_atten;
endif
endfunction
