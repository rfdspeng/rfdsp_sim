% Author: Ryan Tsai
% Created: 2022-02-06

function [c,cinfo] = dpd_training(cfg,x,y)
% cfg.kernel_mask: Nx3 array, where N is the number of LUTs
  % Column 1: Signal delay
  % Column 2: Envelope delay
  % Column 3: 8-bit number. Indicates which envelope orders are active. 0-7 for aligned terms, 1-8 for leading/lagging terms.
% x = PA input, scaled so max(|x|) = 1
% y = PA output, scaled so max(|y|) = 1

if isrow(x), x = x.'; endif
if isrow(y), y = y.'; endif

[Y,Yinfo] = generate_kernel_matrix(cfg.kernel_mask,y);
U = chol(Y'*Y);
L = U';
c = Y\x;
cinfo = Yinfo;

endfunction

function [Y,Yinfo] = generate_kernel_matrix(kernel_mask,y)
env = abs(y);
Y = [];
Yinfo = []; % each row corresponds to a kernel. Indicates signal delay, envelope delay, and envelope order.
for ldx = 1:size(kernel_mask,1)
  lut = kernel_mask(ldx,:);
  delay = lut(1);
  if delay >= 0, yl = [zeros(delay,1); y(1:end-delay)];
  else, yl = [x(1+delay:end); zeros(delay,1)]; endif
  
  env_delay = lut(2);
  if env_delay >= 0, envl = [zeros(env_delay,1); env(1:end-env_delay)];
  else, envl = [env(1+env_delay:end); zeros(env_delay,1)]; endif
  
  mask = dec2bin(lut(3)); mask = fliplr(mask);
  for mdx = 1:length(mask)
    b = mask(mdx);
    if strcmp(b,'0'), continue; endif
    if delay == env_delay, env_order = mdx-1;
    else, env_order = mdx; endif
    
    basis = yl.*envl.^env_order;  
    Y = [Y basis];
    Yinfo = [Yinfo; [delay env_delay env_order]];
  endfor
endfor

endfunction
