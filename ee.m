function [ out ] = ee( I1,I2,n,m )
h=zeros(n,m);
for ii=1:n
    for jj=1:m
        if (I1==ii && I2==jj)
            h(ii,jj)=1;
        else
            h(ii,jj)=0;
        end
    end
end
out=h;
end
