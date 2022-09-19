% Author: Ryan Tsai
% Created: 2022-03-05

function noise_dbc = calculate_noise_dbc(p,f,sl,sh,nl,nh,en_plot)
% Calculate noise relative to signal in dBc
% p = psd in V^2
% f = psd frequencies
% sl, sh = signal low, signal high
% nl, nh = noise low, noise high

if ~exist('en_plot','var'), en_plot = 0; endif

sdx = f >= sl & f <= sh;
ndx = f >= nl & f <= nh;
psig = p(sdx); pnoise = p(ndx);
noise_dbc = 10*log10(sum(pnoise)/sum(psig));

if en_plot
  fsig = f(sdx); fnoise = f(ndx);
  figure; hold on;
  plot(f,10*log10(p),'LineWidth',10);
  plot(fsig,10*log10(psig),'LineWidth',5);
  plot(fnoise,10*log10(pnoise),'LineWidth',5);
  legend("Full PSD","Signal","Noise");
endif

endfunction
