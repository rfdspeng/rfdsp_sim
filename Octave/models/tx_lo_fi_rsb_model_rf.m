% Author: Ryan Tsai
% Created: 2022-04-14

function y = tx_lo_fi_rsb_model_rf(cfg,x)
% Tx FI RSB introduced in LO signals
% This is the RF model
% cfg.gi/gq: I/Q path gain
% cfg.phii/phiq: I/Q phase in degrees
% cfg.fc: LO frequency (MHz)
% cfg.fs: Sampling rate (MHz)

gi = cfg.gi; gq = cfg.gq;
phii = cfg.phii; phiq = cfg.phiq;
fc = cfg.fc; fs = cfg.fs;

wc = fc*2*pi/fs; % digital LO frequency
phii = phii*pi/180; phiq = phiq*pi/180; % convert to radians

% I and Q
len = length(x);
i = real(x); q = imag(x);

% LOs
loi = gi*cos(wc*(0:len-1)+phii);
loq = gq*sin(wc*(0:len-1)+phiq);

y = loi.*i - loq.*q; % RF output
endfunction
