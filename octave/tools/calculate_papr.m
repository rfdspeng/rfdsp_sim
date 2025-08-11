% Author: Ryan Tsai
% Created: 2022-04-14

function [papr,avg_power,peak_power] = calculate_papr(x,p)
% Calculates the p-percentile PAPR of x
% p is in %, e.g. p = 99.99

env2 = x.*conj(x);
env2 = sort(env2);
avg_power = mean(env2);
peak_power = env2(ceil(length(env2)*p/100));
papr = 10*log10(peak_power/avg_power);

endfunction
