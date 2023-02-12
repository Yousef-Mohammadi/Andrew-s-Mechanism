function [ in_m ] = Modify( in )
nn=strfind(in,'--');
np=strfind(in,'-+');
pn=strfind(in,'+-');
% -+ (neg plus):
in_m=[];
for ii=1:length(nn)
    if ii==1 
        in_m=[in_m,in(1:nn(ii)-1),'+'];
        if length(nn)==1
            in_m=[in_m,in(nn(ii)+2:end)];
        end
    elseif ii~=length(nn)
        in_m=[in_m,in(nn(ii-1)+2:nn(ii)-1),'+'];
    elseif ii==length(nn)
        in_m=[in_m,in(nn(ii-1)+2:nn(ii)-1),'+'];
        in_m=[in_m,in(nn(ii)+2:end)];
    end
end
if ~isempty(in_m)
    in=in_m;
    in_m=[];
end
% +- (neg plus):
for ii=1:length(np)
    if ii==1 
        in_m=[in_m,in(1:np(ii)-1),'-'];
        if length(np)==1
            in_m=[in_m,in(np(ii)+2:end)];
        end
    elseif ii~=length(np)
        in_m=[in_m,in(np(ii-1)+2:np(ii)-1),'-'];
    elseif ii==length(np)
        in_m=[in_m,in(np(ii-1)+2:np(ii)-1),'-'];
        in_m=[in_m,in(np(ii)+2:end)];
    end
end
if ~isempty(in_m)
    in=in_m;
    in_m=[];
end
% +- (plus neg):
for ii=1:length(pn)
    if ii==1
        in_m=[in_m,in(1:pn(ii)-1),'-'];
        if length(pn)==1
            in_m=[in_m,in(pn(ii)+2:end)];
        end
    elseif ii~=length(pn)
        in_m=[in_m,in(pn(ii-1)+2:pn(ii)-1),'-'];
    else
        in_m=[in_m,in(pn(ii-1)+2:pn(ii)-1),'-'];
        in_m=[in_m,in(pn(ii)+2:end)];
    end
end
if isempty(in_m)
    in_m=in;
end
if in_m(1)=='+'
    in_m=in_m(2:end);
end