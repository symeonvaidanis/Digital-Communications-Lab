%% Bandpass FIR filter – closed-form and windowing
clear all; close all;
load sima; % Fs=8192;
f1=800; f2=1600; 
Ts=1/Fs;
f2m1=(f2-f1); f2p1=(f2+f1)/2; N=256; %
t=[-(N-1):2:N-1]*Ts/2;
hbp=2/Fs*cos(2*pi*f2p1*t).*sin(pi*f2m1*t)/pi./t;
hbpw=hbp.*kaiser(length(hbp),5)';
wvtool(hbp,hbpw);
pause
sima_bp=conv(s,hbp);
figure; pwelch(sima_bp,[],[],[],Fs);
pause
sima_bpw=conv(s,hbpw);
figure; pwelch(sima_bpw,[],[],[],Fs);
pause