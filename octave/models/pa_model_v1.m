% Author: Ryan Tsai
% Created: 2022-02-06

function y = pa_model(x)
% Simple wrapper function for rapp_saleh_model()
% Parameters specified so that when |x| = 1, the peak value, compression is roughly 2.5dB

% Rapp params
cfg.g = 10^(30/20);
cfg.smoothness = 2;
cfg.osat = 25;

% Saleh params
cfg.a = 0.5; cfg.b = 10;

y = rapp_saleh_model_v1(cfg,x);

endfunction
