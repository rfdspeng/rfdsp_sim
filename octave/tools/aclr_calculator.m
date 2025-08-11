% Author: Ryan Tsai
% Created: 2022-05-06

function [aclr,aclrm,aclrp] = aclr_calculator(cfg,x)

bw = cfg.bw; scs = cfg.scs; fs = cfg.fs;
start_sc = cfg.start_sc; num_sc = cfg.num_sc;

% Integration limits
nrb = bw*5; maxsc = nrb*12;
sigl = -nrb*12*scs/1000/2 + start_sc*scs/1000;
sigh = sigl + (num_sc-1)*scs/1000;
mbw = maxsc*scs/1000;

% PSD
[p,f] = calculate_psd(x,fs,scs/1000);

aclrm = -calculate_noise_dbc(p,f,sigl,sigh,-bw-mbw/2,-bw+mbw/2);
aclrp = -calculate_noise_dbc(p,f,sigl,sigh,+bw-mbw/2,+bw+mbw/2);
aclr = min(aclrm,aclrp);

endfunction
