%% Bandpass FIR filter – Equiripple design (Parks Mc Clellan)
clear all; close all;
load sima; % Fs=8192;
f1=800; f2=1600;
f=2*[0 f1*0.95 f1*1.05 f2*0.95 f2*1.05 Fs/2]/Fs;
hbp_pm=firpm(128, f, [0 0 1 1 0 0]);
% figure; freqz(hpm,1);
wvtool(hbp_pm);
sima_bp=conv(s,hbp_pm);
figure; pwelch(sima_bp,[],[],[],Fs);
