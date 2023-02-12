function [  ] = Converter_str2num( in )
n_for=size(in);
len=length(n_for);
fid=fopen('E:\Andrew system\simple_str2num.m','w');
fprintf(fid,sprintf('function [out]=simple_str2num(in)\n'));
fprintf(fid,sprintf('in=cell(in);\n'));
fprintf(fid,sprintf('out=in;\n'));
I=sprintf('i%d,',[1:len]);
I=I(1:end-1);
for ii=1:len
    fprintf(fid,sprintf('for i%d=1:%d\n',ii,n_for(ii)));
end
fprintf(fid,sprintf('if (~isempty(in{%s}))\n',I));
fprintf(fid,sprintf('if (~isnan(str2double(in{%s})))\n',I));
fprintf(fid,sprintf('out{%s}=str2double(in{%s});\n',I,I));
fprintf(fid,'end\n');
fprintf(fid,'end\n');
for ii=1:len
    fprintf(fid,sprintf('end\n'));
end
fclose(fid);