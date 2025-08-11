% Author: Ryan Tsai
% Created: 2022-03-26

function x = mod_ask(nsym,bps,osr)
% Generates a vector of ASK symbols
% nsym = number of symbols
% bps = bits per symbol (any positive integer)
% osr = oversampling ratio

bitstream = (rand(1,nsym*bps)-0.5) > 0;
valid_syms = -2^bps+1:2:2^bps-1;

x = zeros(1,nsym); xdx = 1;
for bdx = 1:bps:length(bitstream)
  bits = bitstream(bdx:bdx+bps-1);
  idx = 1;
  for bbdx = 1:length(bits)
    idx = idx+bits(bbdx)*2^(bbdx-1);
  endfor
  sym = valid_syms(idx);
  x(xdx) = sym;
  xdx = xdx+1;
endfor

x = upsample(x,osr);
x = filter(ones(1,osr),1,x);
%{
hannwin = hanning(osr).';
for xdx = 1:osr:length(x)
  x(xdx:xdx+osr-1) = x(xdx:xdx+osr-1).*hannwin;
endfor
%}

endfunction
