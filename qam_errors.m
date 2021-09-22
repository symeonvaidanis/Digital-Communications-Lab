function errors=qam_errors(w,Nsymb,nsamp,EbNo,rolloff,T)
% Εξομοιώνεται το πλήρες σύστημα πομπού-δέκτη QAM, με αφετηρία
% δυαδική ακολουθία προς εκπομπή. Γίνεται κωδικοποίηση
% Gray, ενώ Γκαουσιανός θόρυβος προστίθεται στο ζωνοπερατό
% QAM σήμα για δοσμένο σηματοθορυβικό λόγο Eb/No.
% Δίνεται δυνατότητα επιλογής φίλτρου μορφοποίησης βασικής ζώνης
% (ορθογωνικό ή Nyquist), ενώ μετράται το BER μετά τη βέλτιστη
% αποδιαμόρφωση με το αντίστοιχο προσαρμοσμένο φίλτρο.
% M (=2^k, k άρτιο)είναι το μέγεθος του σηματικού αστερισμού,
% mapping είναι το μιγαδικό διάνυσμα των Μ QAM συμβόλων,
% σε διάταξη κωδικοποίησης Gray
% Nsymb είναι το μήκος της ακολουθίας QAM,
% nsamp είναι ο συντελεστής υπερδειγμάτισης, (# δειγμάτων/Τ)
% EbNo είναι ο λόγος Eb/No, in db
%
L=sqrt(2^w); % LxL σηματικός αστερισμός, L=2^l, l>0
l=log2(L); k=2*l; M=L^2;
fc=k; % συχνότητα φέροντος, πολλαπλάσιο του Baud Rate (1/T)
SNR=EbNo-10*log10(nsamp/k/2); % SNR ανά δείγμα σήματος
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
core=[1+i;1-i;-1+i;-1-i];
mapping=core;
if(l>1)
    for j=1:l-1
        mapping=mapping+j*2*core(1);
        mapping=[mapping;conj(mapping)];
        mapping=[mapping;-conj(mapping)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ΠΟΜΠΟΣ %%%%%%%%%%%%
x=floor(2*rand(k*Nsymb,1)); % τυχαία δυαδική ακολουθία
xsym=bi2de(reshape(x,k,length(x)/k).','left-msb')';
y=[];
for n=1:length(xsym)
    y=[y mapping(xsym(n)+1)];
end
%% Ορισμός φίλτρου μορφοποίησης Nyquist
delay = T; % Group delay (# περιόδων Τ)
filtorder = delay*nsamp*2;
shaping_filter = rcosine(1,nsamp,'fir/sqrt',rolloff,delay);
%% Εκπεμπόμενο σήμα
ytx=upsample(y,nsamp);
ytx = conv(ytx,shaping_filter);
% quadrature modulation
m=(1:length(ytx));
s=real(ytx.*exp(1j*2*pi*fc*m/nsamp));
% ---------------------
% Προσθήκη λευκού γκαουσιανού θορύβου
snoisy=awgn(s,SNR,'measured'); % θορυβώδες σήμα
clear ytx xsym s n; % για εξοικονόμηση μνήμης
%%%%%%%%% ΔΕΚΤΗΣ %%%%%%%%%%%%%
% Αποδιαμόρφωση
yrx=2*snoisy.*exp(-1j*2*pi*fc*m/nsamp); clear s;
yrx = conv(yrx,shaping_filter);
yrx = downsample(yrx,nsamp); % Υποδειγμάτιση στο πλέγμα nT.
yrx = yrx(2*delay+(1:length(y))); % περικοπή άκρων συνέλιξης.
% ----------------------
yi=real(yrx); yq=imag(yrx); % συμφασική και εγκάρσια συνιστώσα
xrx=[]; % διάνυσμα δυαδικής ακολουθίας εξόδου -- αρχικά κενό
q=[-L+1:2:L-1];
for n=1:length(yrx) % επιλογή πλησιέστερου σημείου
    [m,j]=min(abs(q-yi(n)));
    yi(n)=q(j);
    [m,j]=min(abs(q-yq(n)));
    yq(n)=q(j);
    m=1;
    while (mapping(m)~=yi(n)+i*yq(n)) m=m+1; end
    xrx=[xrx; de2bi(m-1,k,'left-msb')'];
end
%% ΕΚΤΙΜΗΣΗ ΛΑΘΩΝ
% Με δύο τρόπους: (α) με σύγκριση των σημείων QAM (y-yrx)
% (β) με σύγκριση των δυαδικών ακολουθιών (x-xrx)
errors=sum(not(xrx==x));
end