function [ TT ] = sign1( in )
TT=zeros(1,length(in));
for ii=1:length(in)
    if (in(ii)>0)
        TT(ii)=1;
    else
        TT(ii)=0;
    end
end