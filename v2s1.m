function [ out ] = v2s1( in )
out=0;
in=double(in);
for ii=1:length(in)
    if (in(ii)==1)
        out=1;
        break;
    end
end
end