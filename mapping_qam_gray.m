function mapping_qam_gray(L)
k=log2(L);
core=[1+1j 1-1j -1+1j -1-1j];
mapping=core;
if(k>1)
    for m=1:(k/2-1)
        mapping=mapping+m*1j*2*core(2);
        mapping=[mapping conj(mapping)];
        mapping=[mapping -conj(mapping)];
    end
end
scatterplot(mapping);
for i=1:L
    text(real(mapping(1,i)), imag(mapping(1,i)), num2str(de2bi(i-1,k,'left-msb')), 'FontSize', 10);
end
end