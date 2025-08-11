% Author: Ryan Tsai
% Created: 2022-02-06

function y = rapp_saleh_model(cfg,x)
% https://www.mathworks.com/help/comm/ref/memorylessnonlinearity.html
% Rapp model for AMAM
% Saleh model for AMPM

% Rapp parameters
% cfg.g: Gain in the linear region
% cfg.smoothness: Smoothness factor
% cfg.osat: Output saturation level

% Saleh parameters
% cfg.a: AMPM alpha
% cfg.b: AMPM beta

g = cfg.g; s = 2*cfg.smoothness; osat = cfg.osat;
a = cfg.a; b = cfg.b;

env = abs(x);
ph = angle(x);
env_y = g*env./(1 + (g*env/osat).^s).^(1/s);
ampm = a*env.^2./(1 + b*env.^2); %ampm = 0;
y = env_y.*exp(1j*ph).*exp(1j*ampm);

endfunction
