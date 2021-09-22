function mapping_psk_gray(L)
k=log2(L);
ph1=[pi/4];
theta=[ph1 -ph1 pi-ph1 -pi+ph1];
mapping=exp(1j*theta);
if(k>2)
    for m=3:k
        theta=theta/2;
        help1=mapping.*exp(-1j*theta);
        help2=help1-2*real(help1);
        mapping=[help1 help2];
        theta=angle(mapping);
    end
end
scatterplot(mapping);
for i=1:L
    text(real(mapping(1,i)), imag(mapping(1,i)), num2str(de2bi(i-1,k,'left-msb')), 'FontSize', 10);
end
end