function [ out ] = FU( n,m )
S=zeros(m*n,m*n);
for ii=1:n
    for jj=1:m
        S=S+Diadic(ee(ii,jj,n,m),ee(jj,ii,m,n));
    end
end
out=S;
end

