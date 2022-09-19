% Is halfband better than nonhalfband?

clear; clc; close all;

fs = 30.72;
obw = fs/2;
ripple_1x_spec = 0.01; % dB
rej_1x_spec = 80; % dB

norm_bw = obw/fs;
wp = norm_bw/2;
ws = 1-wp;
f = [0 wp ws 1]; a = [1 1 0 0];
weight = [1 1];

%devp = min([10^(ripple_1x_spec/20) 10^(-ripple_1x_spec/20)]-1);
%devs = 10^(-rej_1x_spec/20);
%dev = [devp devs];
%weight = max(dev)./dev;
%weight = weight.*weight;

n = 2;
while 1
  b = firls(n,f,a,weight); b = b/sum(b);
  [h,w] = freqz(b,1,2^17); w = w/pi;
  hp = h(w <= wp); hs = h(w >= ws);
  ripple_1x = max(abs(20*log10(abs(hp))));
  rej_1x = min(-20*log10(abs(hs)));
  
  passflag = 0;
  if ripple_1x <= ripple_1x_spec/2 && rej_1x >= rej_1x_spec, passflag = 1; endif
  
  if passflag, break; endif
  n = n+2; 
endwhile

rs = 6;
while 1
  b_fx = round(b*2^rs);
  b1 = b_fx/2^rs; b1 = b1/sum(b1);
  [h,w] = freqz(b1,1,2^17); w = w/pi;
  hp = h(w <= wp); hs = h(w >= ws);
  ripple_1x = max(abs(20*log10(abs(hp))));
  rej_1x = min(-20*log10(abs(hs)));
  
  passflag = 0;
  if ripple_1x <= ripple_1x_spec/2 && rej_1x >= rej_1x_spec, passflag = 1; endif
  
  if passflag, break; endif
  rs = rs+1;
endwhile
ntaps = length(unique(b_fx(b_fx~=0)));