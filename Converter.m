function [  ] = Converter( in )
n_for=size(in);
len=length(n_for);
fid=fopen('E:\Andrew system\simple.m','w');
fprintf(fid,sprintf('function [out]=simple(in)\n'));
fprintf(fid,sprintf('in=cell(in);\n'));
I=sprintf('i%d,',[1:len]);
I=I(1:end-1);
for ii=1:len
    fprintf(fid,sprintf('for i%d=1:%d\n',ii,n_for(ii)));
end
fprintf(fid,sprintf('out{%s}=char(vpa(simplify(sym(in{%s})),4));\n',I,I));
for ii=1:len
    fprintf(fid,sprintf('end\n'));
end
fclose(fid);