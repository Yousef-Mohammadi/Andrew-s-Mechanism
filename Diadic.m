function [ C ] = Diadic( A,B )
ar=size(A,1); ac=size(A,2);
br=size(B,1); bc=size(B,2);
for ii=1:ar
    for jj=1:ac
        C((ii-1)*br+1:ii*br,(jj-1)*bc+1:jj*bc)=A(ii,jj)*B;
    end
end
end

