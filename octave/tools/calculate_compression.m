% Author: Ryan Tsai
% Created: 2022-05-03

function comp = calculate_compression(x,y,en_plot,titlein)
% x and y are input and output signals
% comp is compression in dB

if ~exist('en_plot','var'), en_plot = 0; endif

if isrow(x), x = x.'; endif
if isrow(y), y = y.'; endif
x = x-mean(x); y = y-mean(y);
x = abs(x)/max(abs(x)); y = abs(y)/max(abs(y));

ks = 1:9; X = [];
for k = ks
  X = [X x.^k];
endfor
c = X\y;

ye = 0; xe = linspace(0,1,2^16);
for kdx = 1:length(ks);
  k = ks(kdx);
  ye = ye + c(kdx)*xe.^k;
endfor

g = 20*log10(ye(2:end)./xe(2:end)); g = g-max(g);
comp = -g(end);

if en_plot
  figure; hold on;
  plot(x,y,'.','markersize',20);
  plot(xe,ye,'linewidth',15);
  if exist('titlein','var'), title(titlein);
  else, title("AMAM"); endif
  xlabel("Normalized Input Envelope"); ylabel("Normalized Output Envelope");
  grid on;
  legend("Data","Polyfit",'location','south');
  set(gca,'fontsize',40);

  figure;
  plot(xe(2:end),g,'.','markersize',20);
  if exist('titlein','var'), title(titlein);
  else, title("Normalized Gain"); endif
  xlabel("Normalized Input Envelope"); ylabel("Normalized Gain (dB)");
  grid on;
  set(gca,'fontsize',40);
endif

endfunction
