% Author: Ryan Tsai
% Created: 2022-04-14

function y = tx_lo_fi_rsb_model_bb(cfg,x)
% Tx FI RSB introduced in LO signals
% This is the BB model
% cfg.gi/gq: I/Q path gain
% cfg.phii/phiq: I/Q phase in degrees

gi = cfg.gi; gq = cfg.gq;
phii = cfg.phii; phiq = cfg.phiq;
phii = phii*pi/180; phiq = phiq*pi/180; % convert to radians

% I and Q
len = length(x);
i = real(x); q = imag(x);

yi = gi*cos(phii)*i - gq*sin(phiq)*q;
yq = gi*sin(phii)*i + gq*cos(phiq)*q;

y = yi + 1j*yq;
endfunction
