% Author: Ryan Tsai
% Created: 2022-05-03

function y = cfr(x,target_papr)
% CFR model (just hard clipping for now)

v = rms(x); % rms voltage
th = v*10^(target_papr/20);
env = abs(x);
ph = angle(x);
env(env > th) = th;
y = env.*exp(1j*ph);

endfunction
