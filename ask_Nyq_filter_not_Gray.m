function errors=ask_Nyq_filter_not_Gray(k,Nsymb,nsamp,EbNo,rolloff,T)
%% ΠΑΡΑΜΕΤΡΟΙ 
% k είναι ο αριθμός των bits ανά σύμβολο, έτσι L=2^k
% Nsymb είναι το μήκος της εξομοιούμενης ακολουθίας συμβόλων L-ASK
% nsamp είναι ο συντελεστής υπερδειγμάτισης, δηλ. #samples/Td
% EbNo είναι ο ανηγμένος σηματοθορυβικός λόγος, Εb/No, σε db
% Συντελεστής πτώσης -- rolloff factor
% T τέτοιο ώστε group delay = T*nsamp
L=2^k; step=2; mapping=[-L+1:step:L-1];
SNR=EbNo-10*log10(nsamp/2/k); % SNR ανά δείγμα σήματος
z=randi([0 1], 1,Nsymb*k);
zsym=bi2de(reshape(z,k,length(z)/k).','left-msb');
x=[]; 
for i=1:length(zsym) 
    x=[x mapping(zsym(i)+1)]; 
end
%% Ορισμός παραμέτρων φίλτρου
delay=T; % Group delay (# of input symbols)
filtorder=delay*nsamp*2; % τάξη φίλτρου
% κρουστική απόκριση φίλτρου τετρ. ρίζας ανυψ. συνημιτόνου
rNyquist=rcosine(1,nsamp,'fir/sqrt',rolloff,delay);
% ----------------------
%% ΕΚΠΕΜΠΟΜΕΝΟ ΣΗΜΑ
% Υπερδειγμάτιση και εφαρμογή φίλτρου rNyquist
y=upsample(x,nsamp);
ytx = conv(y,rNyquist);
ynoisy=awgn(ytx,SNR,'measured'); % θορυβώδες σήμα
% ----------------------
%% ΛΑΜΒΑΝΟΜΕΝΟ ΣΗΜΑ
% Φιλτράρισμα σήματος με φίλτρο τετρ. ρίζας ανυψ. συνημ.
yrx=conv(ynoisy,rNyquist);
yrx = downsample(yrx,nsamp); % Υποδειγμάτιση
yrx = yrx(2*delay+1:end-2*delay); % περικοπή, λόγω καθυστέρησης
% Ανιχνευτής ελάχιστης απόστασης L πλατών
w_bits=[];
for i=1:length(yrx)
    [m,j]=min(abs(mapping-yrx(i)));
    w_bits=[w_bits de2bi(j-1,k,'left-msb')];
end
err=not(z==w_bits);
errors=sum(err);
end