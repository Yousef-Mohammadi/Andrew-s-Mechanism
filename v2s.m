function [ out ] = v2s( in )
out=1;
in=double(in);
for ii=1:length(in)
    if (in(ii)==0)
        out=0;
        break;
    end
end
end

