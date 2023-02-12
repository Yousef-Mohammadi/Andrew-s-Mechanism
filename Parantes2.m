function out  = Parantes2( in,tag )
plus=strfind(in,'+');
neg=strfind(in,'-');
mul=strfind(in,'*');
div=strfind(in,'/');
pari=strfind(in,'(');
parf=strfind(in,')');
Pari=[]; Parf=[];
for ii=1:length(pari)
    if (length(find(pari<pari(ii)))==length(find(parf<pari(ii))))
        Pari=[Pari,pari(ii)];
    end
end
for ii=1:length(parf)
    if (length(find(pari<parf(ii)))==length(find(parf<=parf(ii))))
        Parf=[Parf,parf(ii)];
    end
end
Plus=[]; Neg=[]; Mul=[]; Div=[];
for ii=1:length(plus)
    cc=0;
    for jj=1:length(Pari)
        if (plus(ii)<Pari(jj) || plus(ii)>Parf(jj))
            cc=cc+1;
        end
    end
    if cc==length(Pari)
        Plus=[Plus,plus(ii)];
    end
end
for ll=1:length(Plus)
    if Plus(ll)>1
        if in(Plus(ll)-1)=='e'
            Plus(ll)=nan;
        end
    end
end
Plus(isnan(Plus))=[];
for ii=1:length(neg)
    cc=0;
    for jj=1:length(Pari)
        if (neg(ii)<Pari(jj) || neg(ii)>Parf(jj))
            cc=cc+1;
        end
    end
    if cc==length(Pari)
        Neg=[Neg,neg(ii)];
    end
end
for ll=1:length(Neg)
    if Neg(ll)>1
        if in(Neg(ll)-1)=='e'
            Neg(ll)=nan;
        end
    end
end
Neg(isnan(Neg))=[];
for ii=1:length(mul)
    cc=0;
    for jj=1:length(Pari)
        if (mul(ii)<Pari(jj) || mul(ii)>Parf(jj))
            cc=cc+1;
        end
    end
    if cc==length(Pari)
        Mul=[Mul,mul(ii)];
    end
end
for ii=1:length(div)
    cc=0;
    for jj=1:length(Pari)
        if (div(ii)<Pari(jj) || div(ii)>Parf(jj))
            cc=cc+1;
        end
    end
    if cc==length(Pari)
        Div=[Div,div(ii)];
    end
end
if (~strcmp(tag,'n') && ~strcmp(tag,'^'))
    if ~isempty(Neg)
        if Neg(1)==1
            Neg=Neg(2:end);
        end
    end
end
if ~isempty(Plus)
    if Plus(1)==1
        Plus=Plus(2:end);
    end
end
if (~isempty(Neg) || ~isempty(Plus) || ~isempty(Mul) || ~isempty(Div))
    out=sprintf('(%s)',in);
else
    out=in;
end