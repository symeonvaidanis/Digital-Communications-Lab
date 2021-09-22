%% Lowpass FIR filter – closed-form and windowing
clear all; close all;
load sima; % Fs=8192;
fo=1200; Ts=1/Fs;
N=128; % Τάξη φίλτρου)
t=[-(N-1):2:N-1]*Ts/2;
hlp=1/Fs*sin(2*pi*fo*t)/pi./t;
hlpw=hlp.*kaiser(length(hlp),5)';
wvtool(hlp,hlpw);
pause
sima_lp=conv(s,hlpw);
figure; pwelch(sima_lp,[],[],[],Fs);
pause
sima_lp=conv(s,hlp);
figure; pwelch(sima_lp,[],[],[],Fs);
pause