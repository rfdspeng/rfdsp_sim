% Author: Ryan Tsai
% Created: 2022-04-02

% If you assume all symbols are the same length, then you can also vectorize the symbol slicing
% And if you assume no dynamic power change, you can also vectorize the single-tap TD equalizer

% The standard FD equalizer and my equalizer are the same - how?

function evm = ofdm_evm_calculator(cfg,x,y)
% Calculates EVM (%) for OFDM waveform
% x = reference, y = signal
% Signals must be at baseband sampling rate and time-aligned

if ~isfield(cfg,'en_plot'), en_plot = 0;
else, en_plot = cfg.en_plot; endif
if ~isfield(cfg,'en_fd_eq'), en_fd_eq = 0;
else, en_fd_eq = cfg.en_fd_eq; endif

% Convert from time domain symbols to frequency domain subcarriers
n = min(length(x),length(y));
nfft = cfg.nfft; ncp = cfg.ncp; sym_len = nfft+ncp;
fft_start = max(round(ncp/2),1);
tone_idx = cfg.tone_idx;
tone_idx = tone_idx == 2;
tones_xx = []; % nsym x ntone
tones_yy = []; % nsym x ntone
for sdx = 1:sym_len:n-sym_len+1
  sym_x = x(sdx:sdx+sym_len-1);
  sym_y = y(sdx:sdx+sym_len-1);

  % Single-tap time-domain equalizer
  yy = sym_y.'; xx = sym_x.';
  xx_hat = (xx'*yy)/(yy'*yy)*yy;
  sym_y = xx_hat.';

  % Slice symbols and circshift
  sym_x = sym_x(fft_start:fft_start+nfft-1);
  sym_y = sym_y(fft_start:fft_start+nfft-1);
  nshift = ncp-fft_start+1;
  sym_x = circshift(sym_x,-nshift);
  sym_y = circshift(sym_y,-nshift);

  % FFT
  tones_x = fft(sym_x); tones_y = fft(sym_y);
  tones_x = tones_x(tone_idx);
  tones_y = tones_y(tone_idx);
  tones_xx = [tones_xx; tones_x];
  tones_yy = [tones_yy; tones_y];
endfor

% Frequency domain equalization
fd_eqs = ones(size(tones_x));
if en_fd_eq
  for edx = 1:length(fd_eqs)
    xx = tones_xx(:,edx); yy = tones_yy(:,edx);

    % Standard definition

    num = xx'*xx; den = xx'*yy; % NS'*NS and NS'*MS
    fd_eq = num/den; % (NS'*NS)/(NS'*MS)
    yyy = fd_eq*yy; % equalize
    %}

    % My definition - this also works
    %{
    %num = xx'*yy; den = yy'*yy; % does not work -> dot product of complex vectors is not commutative
    num = yy'*xx; den = yy'*yy;
    fd_eq = num/den;
    yyy = fd_eq*yy;
    %}

    %{
    figure; hold on;
    plot(real(xx),imag(xx),'x');
    plot(real(yy),imag(yy),'x');
    plot(real(yyy),imag(yyy),'x');
    %}

    tones_yy(:,edx) = yyy;
    fd_eqs(edx) = fd_eq;
  endfor
endif

% Inverse transform precoding
if cfg.en_tprecode
  for sdx = 1:size(tones_xx,1)
    tones_xx(sdx,:) = ifft(tones_xx(sdx,:));
    tones_yy(sdx,:) = ifft(tones_yy(sdx,:));
  endfor
endif

% Calculate EVM
err = tones_xx-tones_yy;
errv = reshape(err.',1,numel(err));
tonesv = reshape(tones_xx.',1,numel(tones_xx));
errp = errv*errv';
tonesp = tonesv*tonesv';
evm = sqrt(errp/tonesp)*100;

if en_plot
  if isfield(cfg,'title'), titlestr = cfg.title;
  else, titlestr = "Constellation"; endif
  tones_x = reshape(tones_xx.',1,numel(tones_xx));
  tones_y = reshape(tones_yy.',1,numel(tones_yy));
  a_norm = rms(tones_x);
  tones_x = tones_x/a_norm;
  tones_y = tones_y/a_norm;
  figure; hold on;
  plot(tones_x,'x','markersize',25);
  plot(tones_y,'x','markersize',10);
  title(titlestr);
  xlabel("I"); ylabel("Q");
  set(gca,"fontsize",25);
endif

endfunction
