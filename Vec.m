function [ out ] = Vec( A )
n=size(A,1);
m=size(A,2);
for ii=1:m
   out((ii-1)*n+1:ii*n,1)=A(:,ii);
end
end

