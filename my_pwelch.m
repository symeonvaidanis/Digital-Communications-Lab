function my_pwelch(x,Fs)
%% PART A (One side Power Spectral Density)
x=x';
A=size(x);
L=A(1,2);
N=2^nextpow2(L);
Fo=Fs/N;
f=(0:N-1)*Fo;
X=fft(x,N);
power=X.*conj(X)/N/L;
figure(1)
plot(f(1:(N/2)),power(1:(N/2)))
xlabel('Frequency (Hz)') 
ylabel('Power') 
title('bf The oneside power of x - VAIDANIS 18005') 
%% PART B (Find how many parts of signal i will take)
K=2^nextpow2(L/8);
if K < 256
    K = 256;
end
G=0; %the number of parts of thw signal
for i=K:(K/2):L
    G=G+1;
end
s=zeros(G,K); %each row is a part
for i=1:G
    k=(i-1)*K/2;
    for j=1:K
        s(i,j)=x(1,j+k);
    end
end
%% PART C (FFT for each row, periodiagramm for eack row)
S=zeros(G,K); %each column is for a frequency
for i=1:G
    A=fft(s(i,:),K);
    S(i,:)=A.*conj(A)/Fs/K;
end
S_mean=db(mean(S)); %dB and mean value for each frequency
Fonew=Fs/K;
fnew=(0:K-1)*Fonew;
figure(2)
plot(fnew(1:(K/2)),S_mean(1:(K/2)))
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
title('My implementation to estimate PSD of the signal - VAIDANIS 18005')
%% END
disp('Finish')
end