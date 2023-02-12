function [out]=simple(in)
in=cell(in);
for i1=1:1
for i2=1:1
out{i1,i2}=char(vpa(simplify(sym(in{i1,i2})),4));
end
end
