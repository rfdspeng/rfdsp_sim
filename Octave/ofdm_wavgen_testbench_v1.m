% Test ofdm_wavgen()

clear; clc; close all;

%nsym = 140; bw = 20; scs = 15; lcrb = 18; srb = 41; modorder = 4; en_tprecode = 0;
%[x,fs,nfft,ncp,nrb] = ofdm_wavgen(nsym,bw,scs,lcrb,srb,modorder,en_tprecode);

nsym = 140*10; bw = 20; scs = 15; lcsc = 1; ssc = 41*12; modorder = 4; en_tprecode = 0;
[x,fs,nfft,ncp,nrb] = ofdm_wavgen_by_subcarrier(nsym,bw,scs,lcsc,ssc,modorder,en_tprecode);

%[p,f] = pwelch(x,hanning(nfft*16),[],nfft*16,fs);
fs = 2*pi;
[p,f,num_segments] = pwelch_custom(x,hanning(nfft*16),[],nfft*16,fs);
%p = fftshift(p); f = f-fs/2;

%{
sigl = -nrb*12*scs/1000/2 + srb*12*scs/1000;
sigh = sigl + lcrb*12*scs/1000 - scs/1000;
fsig = f(f >= sigl & f <= sigh);
psig = p(f >= sigl & f <= sigh);
figure; hold on;
plot(f,10*log10(p));
plot(fsig,10*log10(psig));
xlim([-fs/2 fs/2]); grid on;
%}

figure; hold on;
plot(f,10*log10(p));
xlim([-fs/2 fs/2]); grid on;

nsym = 140*10; bw = 20; scs = 15; lcsc = 1; ssc = 41*12+1; modorder = 4; en_tprecode = 0;
[x,fs,nfft,ncp,nrb] = ofdm_wavgen_by_subcarrier(nsym,bw,scs,lcsc,ssc,modorder,en_tprecode);
fs = 2*pi;
[p,f,num_segments] = pwelch_custom(x,hanning(nfft*16),[],nfft*16,fs);
plot(f,10*log10(p));

nsym = 140*10; bw = 20; scs = 15; lcsc = 1; ssc = 41*12+2; modorder = 4; en_tprecode = 0;
[x,fs,nfft,ncp,nrb] = ofdm_wavgen_by_subcarrier(nsym,bw,scs,lcsc,ssc,modorder,en_tprecode);
fs = 2*pi;
[p,f,num_segments] = pwelch_custom(x,hanning(nfft*16),[],nfft*16,fs);
plot(f,10*log10(p));

nsym = 140*10; bw = 20; scs = 15; lcsc = 1; ssc = 41*12+3; modorder = 4; en_tprecode = 0;
[x,fs,nfft,ncp,nrb] = ofdm_wavgen_by_subcarrier(nsym,bw,scs,lcsc,ssc,modorder,en_tprecode);
fs = 2*pi;
[p,f,num_segments] = pwelch_custom(x,hanning(nfft*16),[],nfft*16,fs);
plot(f,10*log10(p));

nsym = 140*10; bw = 20; scs = 15; lcsc = 1; ssc = 41*12+4; modorder = 4; en_tprecode = 0;
[x,fs,nfft,ncp,nrb] = ofdm_wavgen_by_subcarrier(nsym,bw,scs,lcsc,ssc,modorder,en_tprecode);
fs = 2*pi;
[p,f,num_segments] = pwelch_custom(x,hanning(nfft*16),[],nfft*16,fs);
plot(f,10*log10(p));