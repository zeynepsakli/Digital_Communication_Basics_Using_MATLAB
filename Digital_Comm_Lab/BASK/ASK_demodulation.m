function [ data_demodule ] = ASK_demodulation( data,carrier,fs,n,Tb)
%UNTÝTLED Summary of this function goes here
%   Detailed explanation goes here
    
    modulated=data.*carrier;
    y=modulated.*carrier;
    KE=(10.^2).*(Tb/4);
    Q=zeros(1,n);
    demodulated=zeros(1,n);
    for i=1:n
        Q(i)=trapz(y(fs*(i-1)+1:fs*i))/fs;
        demodulated(i)=(Q(i)>=KE);
    end
    
    i=1;
    data_demodule=zeros(1,(n*fs)/Tb);
    for j=1:Tb:n
        for i=i:fs+i;
            data_demodule(i)=demodulated(j);
        end
    end
    
end

