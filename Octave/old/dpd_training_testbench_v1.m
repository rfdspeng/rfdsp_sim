clear; clc; close all;

cfg_wav.fs = 491.52; cfg_wav.bw = 100;
cfg_wav.length = 2^17;
x = wavgen(cfg_wav); x = x/max(abs(x));
y = pa_model(x);

fs = cfg_wav.fs; rbw = 15/1000;
[px,f] = calculate_psd(x,fs,rbw);
py = calculate_psd(y,fs,rbw);
figure; hold on;
plot(f,10*log10(px/max(px))); plot(f,10*log10(py/max(py)));
ylim([-100 0]); grid on;

cfg_dpd.kernel_mask = [0 0 255];
                       %1 1 31;
                       %2 2 31
                       %3 3 31];
[c,cinfo] = dpd_training(cfg_dpd,x,y/max(abs(y)));

env = abs(x);
x_dpd = 0;
for cdx = 1:length(c)
  delay = cinfo(cdx,1);
  env_delay = cinfo(cdx,2);
  env_order = cinfo(cdx,3);
  if delay >= 0, xl = [zeros(1,delay) x(1:end-delay)];
  else, xl = [x(1+delay:end) zeros(1,delay)]; endif
  if env_delay >= 0, envl = [zeros(1,env_delay) env(1:end-env_delay)];
  else, envl = [env(1+env_delay:end) zeros(1,env_delay)]; endif
  x_dpd = x_dpd + c(cdx)*xl.*envl.^env_order;
endfor



%x_dpd = c(1)*x + c(2)*x.*env.^2 + c(3)*x.*env.^4;
%x_dpd = c(1)*x + c(2)*x.*env + c(3)*x.*env.^2 + c(4)*x.*env.^3 + c(5)*x.*env.^4;
y_dpd = pa_model(x_dpd);
py_dpd = calculate_psd(y_dpd,fs,rbw);
figure; hold on;
plot(f,10*log10(px/max(px))); plot(f,10*log10(py/max(py))); plot(f,10*log10(py_dpd/max(py_dpd)));
ylim([-100 0]); grid on;