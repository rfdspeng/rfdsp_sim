% Author: Ryan Tsai
% Created: 2022-02-06

function y = wavgen(cfg)
% cfg.length: Number of samples
% cfg.fs: Sampling rate
% cfg.bw: BW

len = cfg.length; fs = cfg.fs; bw = cfg.bw;

x = (2*rand(1,len)-1) + 1j*(2*rand(1,len)-1);
dev = [1e-3 1e-2];
n = kaiserord([bw/fs bw/fs*1.01],[1 0],[1e-3 1e-2]);
b = firls(2*n,[0 bw/fs bw/fs*1.01 1],[1 1 0 0],max(dev)./dev);
x = [x zeros(1,n)];
y = filter(b,1,x);
y = y(1+n:end);

endfunction