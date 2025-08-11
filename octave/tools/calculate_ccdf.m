% Author: Ryan Tsai
% Created: 2022-05-06

function [ccdf,bins] = calculate_ccdf(x,en_plot,titlein)

if ~exist('en_plot','var'), en_plot = 0; endif

bins = -20:0.1:15; % CCDF limits: Instantaneous power to mean power in dB
%bins = -50:0.1:30;
%bins = [-20:0.1:20];

pinst = x.*conj(x);
pavg = mean(pinst);
pratio = 10*log10(pinst/pavg);
histx = histc(pratio,bins);
histx = histx/sum(histx); % normalize to 1
ccdf = 1-cumsum(histx);

if en_plot
  figure;
  semilogy(bins,ccdf);
  xlabel("Inst Power to Mean Power Ratio (dB)");
  ylabel("CCDF (Normalized to 1)");
  if exist('titlein','var'), title(titlein);
  else, title("CCDF"); endif
  set(gca,"fontsize",40);
  xlim([min(bins) max(bins)]);
  ylim([0.0001 1]);
  grid on;
endif

endfunction
