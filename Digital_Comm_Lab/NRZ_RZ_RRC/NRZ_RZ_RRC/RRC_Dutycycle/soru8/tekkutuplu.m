function y=tekkutuplu(isaret,fd,fs,kodlama)

oran=fs/fd;
if strcmp(kodlama,'nrz')
    for i=0:max(size(isaret))-1
        if isaret(i+1)==0
            y(i*oran+1:(i+1)*oran)=zeros(1,oran);
        else
            y(i*oran+1:(i+1)*oran)=ones(1,oran);
        end
    end
elseif strcmp(kodlama,'rz')
    for i=0:max(size(isaret))-1
        if isaret(i+1)==0
            y(i*oran+1:(i+1)*oran)=zeros(1,oran);
        else
            y(i*oran+1:ceil((i+0.5)*oran))=ones(1,ceil(oran*0.5));
            y(ceil((i+0.5)*oran+1):(i+1)*oran)=zeros(1,oran-round(oran*0.5));
        end
    end
else
    error('hatalý hat kodu');
end