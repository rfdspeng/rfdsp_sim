% Check CIC droop for variable R

clear; clc; close all;

norm_bw = 0.5; % normalized signal BW prior to CIC

Rs = 2:64;
%Rs = 2.^(1:6);
%Rs = [2 4];

figure; hold on;
for Rdx = 1:length(Rs)
  R = Rs(Rdx);
  b = ones(1,R)/R;
  [h,w] = freqz(b,1,2^16); w = w/pi;
  hmag = 20*log10(abs(h));
  norm_bw_r = norm_bw/R; % normalized signal BW after zero-stuffing
  idx_sig = w < norm_bw_r;
  w = w(idx_sig)*R; hmag = hmag(idx_sig); % inband response
  plot(w,hmag,'LineWidth',5);
endfor

title("CIC Droop for Variable R");
xlabel("Normalized Digital Frequency");
ylabel("Mag Response (dB)");
set(gca,"FontSize",40);
xlim([0 norm_bw]);
grid on;