% Author: Ryan Tsai
% Created: 2022-03-31

function p = scale_psd(p,f,bw,scs,start_sc,num_sc)
% Scale the PSD so that the average signal bin (in dB) is 0
nrb = bw*5;
sigl = -nrb*12*scs/1000/2 + start_sc*scs/1000;
sigh = sigl + (num_sc-1)*scs/1000;
psig = p(f >= sigl & f <= sigh);
psig = sum(psig)/length(psig);
p = p/psig;
endfunction
