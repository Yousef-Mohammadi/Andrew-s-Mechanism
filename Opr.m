function [ out ] = Opr( a,b,tag )
[m1,n1]=size(a);
[m2,n2]=size(b);
switch tag
    case '+'
        if (m1==1 && n1==1 && m2==1 && n2==1)
            if (isa(a{1},'char') || isa(b{1},'char'))
                if (a{1}==0)
                    out{1,1}=b{1};
                elseif (b{1}==0)
                    out{1,1}=a{1};
                else
                    out{1,1}=Modify(sprintf('%s+%s',num2str(a{1}),num2str(b{1})));
                end
            elseif (isa(a{1},'numeric') && isa(b{1},'numeric'))
                out{1,1}=a{1}+b{1};
            else
                error('Not defined');
            end
        elseif (m1==1 && n1==1)
            for ii=1:m2
                for jj=1:n2
                    if (isa(a{1},'char') || isa(b{ii,jj},'char'))
                        if (a{1}==0)
                            var=b{ii,jj};
                        elseif (b{ii,jj}==0)
                            var=a{1};
                        else
                            var=Modify(sprintf('%s+%s',num2str(a{1}),num2str(b{ii,jj})));
                        end
                    else
                        var=a{1}+b{ii,jj};
                    end
                    out{ii,jj}=var;
                end
            end
        elseif (m2==1 && n2==1)
            for ii=1:m1
                for jj=1:n1
                    if (isa(a{ii,jj},'char') || isa(b{1},'char'))
                        if (a{ii,jj}==0)
                            var=b{1};
                        elseif (b{1}==0)
                            var=a{ii,jj};
                        else
                            var=Modify(sprintf('%s+%s',num2str(a{ii,jj}),num2str(b{1})));
                        end
                    else
                        var=a{ii,jj}+b{1};
                    end
                    out{ii,jj}=var;
                end
            end
        else
            for ii=1:m2
                for jj=1:n2
                    if (m1~=m2 || n1~=n2)
                        error('Matrix dimension do not agree');
                    end
                    if (isa(a{ii,jj},'char') || isa(b{ii,jj},'char'))
                        if (a{ii,jj}==0)
                            var=b{ii,jj};
                        elseif (b{ii,jj}==0)
                            var=a{ii,jj};
                        else
                            var=Modify(sprintf('%s+%s',num2str(a{ii,jj}),num2str(b{ii,jj})));
                        end
                    else
                        var=a{ii,jj}+b{ii,jj};
                    end
                    out{ii,jj}=var;
                end
            end
        end
    case '-'
        if (m1==1 && n1==1 && m2==1 && n2==1)
            if (isa(a{1},'char') || isa(b{1},'char'))
                if (a{1}==0)
                    out{1,1}=Modify(sprintf('-%s',Parantes(b{1})));
                elseif (b{1}==0)
                    out{1,1}=a{1};
                else
                    out{1,1}=Modify(sprintf('%s-%s',num2str(a{1}),Parantes(num2str(b{1}))));
                end
            elseif (isa(a{1},'numeric') && isa(b{1},'numeric'))
                out{1,1}=a{1}-b{1};
            else
                error('Not defined');
            end
        elseif (m1==1 && n1==1)
            for ii=1:m2
                for jj=1:n2
                    if (isa(a{1},'char') || isa(b{ii,jj},'char'))
                        if (a{1}==0)
                            var=Modify(sprintf('-%s',Parantes(b{ii,jj})));
                        elseif (b{ii,jj}==0)
                            var=a{1};
                        else
                            var=Modify(sprintf('%s-%s',num2str(a{1}),Parantes(num2str(b{ii,jj}))));
                        end
                    else
                        var=a{1}-b{ii,jj};
                    end
                    out{ii,jj}=var;
                end
            end
        elseif (m2==1 && n2==1)
            for ii=1:m1
                for jj=1:n1
                    if (isa(a{ii,jj},'char') || isa(b{1},'char'))
                        if (a{ii,jj}==0)
                            var=Modify(sprintf('-%s',Parantes(b{1})));
                        elseif (b{1}==0)
                            var=a{ii,jj};
                        else
                            var=Modify(sprintf('%s-%s',num2str(a{ii,jj}),Parantes(num2str(b{1}))));
                        end
                    else
                        var=a{ii,jj}-b{1};
                    end
                    out{ii,jj}=var;
                end
            end
        else
            for ii=1:m1
                for jj=1:n1
                    if (m1~=m2 || n1~=n2)
                        error('Matrix dimension do not agree');
                    end
                    if (isa(a{ii,jj},'char') || isa(b{ii,jj},'char'))
                        if (a{ii,jj}==0)
                            var=Modify(sprintf('-%s',Parantes(b{ii,jj})));
                        elseif (b{ii,jj}==0)
                            var=a{ii,jj};
                        else
                            var=Modify(sprintf('%s-%s',num2str(a{ii,jj}),Parantes(num2str(b{ii,jj}))));
                        end
                    else
                        var=a{ii,jj}-b{ii,jj};
                    end
                    out{ii,jj}=var;
                end
            end
        end
    case '*'
        if (m1==1 && n1==1 && m2==1 && n2==1)
            if (isa(a{1},'char') || isa(b{1},'char'))
                if (v2s(a{1,1}==0) || v2s(b{1,1}==0))
                    out{1,1}=0;
                elseif (a{1}==1)
                    out{1,1}=b{1};
                elseif (b{1}==1)
                    out{1,1}=a{1};
                elseif (a{1}==-1)
                    out{1,1}=Modify(sprintf('-%s',Parantes(b{1})));
                elseif (b{1}==-1)
                    out{1,1}=Modify(sprintf('-%s',Parantes(a{1})));
                else
                    out{1,1}=sprintf('%s*%s',Parantes(num2str(a{1})),Parantes(num2str(b{1})));
                end
            elseif (isa(a{1},'numeric') && isa(b{1},'numeric'))
                out{1,1}=a{1}*b{1};
            else
                error('Not defined');
            end
        elseif (m1==1 && n1==1)
            for ii=1:m2
                for jj=1:n2
                    if (isa(a{1},'char') || isa(b{ii,jj},'char'))
                        if (v2s(a{1}==0) || v2s(b{ii,jj}==0))
                            out{ii,jj}=0;
                        elseif (a{1}==1)
                            out{ii,jj}=b{ii,jj};
                        elseif (b{ii,jj}==1)
                            out{ii,jj}=a{1};
                        elseif (a{1}==-1)
                            out{ii,jj}=Modify(sprintf('-%s',Parantes(b{ii,jj})));
                        elseif b{ii,jj}==-1
                            out{ii,jj}=Modify(sprintf('-%s',Parantes(a{1})));
                        else
                            out{ii,jj}=sprintf('%s*%s',Parantes(num2str(a{1})),Parantes(num2str(b{ii,jj})));
                        end
                    elseif (isa(a{1},'numeric') && isa(b{ii,jj},'numeric'))
                        out{ii,jj}=a{1}*b{ii,jj};
                    else
                        error('Not defined');
                    end
                end
            end
        elseif (m2==1 && n2==1)
            for ii=1:m1
                for jj=1:n1
                    if (isa(a{ii,jj},'char') || isa(b{1},'char'))
                        if (v2s(a{ii,jj}==0) || v2s(b{1}==0))
                            out{ii,jj}=0;
                        elseif (a{ii,jj}==1)
                            out{ii,jj}=b{1};
                        elseif (b{1}==1)
                            out{ii,jj}=a{ii,jj};
                        elseif a{ii,jj}==-1
                            out{ii,jj}=Modify(sprintf('-%s',Parantes(b{1})));
                        elseif b{1}==-1
                            out{ii,jj}=Modify(sprintf('-%s',Parantes(a{ii,jj})));
                        else
                            out{ii,jj}=sprintf('%s*%s',Parantes(num2str(a{ii,jj})),Parantes(num2str(b{1})));
                        end
                    elseif (isa(a{ii,jj},'numeric') && isa(b{1},'numeric'))
                        out{ii,jj}=a{ii,jj}*b{1};
                    else
                        error('Not defined');
                    end
                end
            end
        else
            if (n1~=m2)
                error('Matrix dimension do not agree');
            end
            for ii=1:m1
                for jj=1:n2
                    num=0; ch=[];
                    for kk=1:n1
                        if ((isa(a{ii,kk},'char') || isa(b{kk,jj},'char')) && v2s1(a{ii,kk}~=0) && v2s1(b{kk,jj}~=0))
                            if (~isempty(ch))
                                if (a{ii,kk}==1)
                                    ch=sprintf('%s+%s',ch,b{kk,jj});
                                elseif (b{kk,jj}==1)
                                    ch=sprintf('%s+%s',ch,a{ii,kk});
                                elseif a{ii,kk}==-1
                                    ch=Modify(sprintf('%s-%s',ch,Parantes(b{kk,jj})));
                                elseif b{kk,jj}==-1
                                    ch=Modify(sprintf('%s-%s',ch,Parantes(a{ii,kk})));
                                else
                                    ch=Modify(sprintf('%s+%s*%s',ch,Parantes(num2str(a{ii,kk})),Parantes(num2str(b{kk,jj}))));
                                end
                            else
                                if (a{ii,kk}==1)
                                    ch=sprintf('%s',num2str(b{kk,jj}));
                                elseif (b{kk,jj}==1)
                                    ch=sprintf('%s',num2str(a{ii,kk}));
                                elseif a{ii,kk}==-1
                                    ch=Modify(sprintf('-%s',Parantes(b{kk,jj})));
                                elseif b{kk,jj}==-1
                                    ch=Modify(sprintf('-%s',Parantes(a{ii,kk})));
                                else
                                    ch=sprintf('%s*%s',Parantes(num2str(a{ii,kk})),Parantes(num2str(b{kk,jj})));
                                end
                            end
                        elseif (isa(a{ii,kk},'numeric') && isa(b{kk,jj},'numeric'))
                            num=num+a{ii,kk}*b{kk,jj};
                        end
                    end
                    if (isempty(ch))
                        out{ii,jj}=num;
                    elseif (num~=0)
                        out{ii,jj}=Modify(sprintf('%d+%s',num,ch));
                    else
                        out{ii,jj}=sprintf('%s',ch);
                    end
                end
            end
        end
    case '-1'
        for ii=1:m1
            for jj=1:m1
                var=a;
                var(ii,:)=[]; var(:,jj)=[];
                var=detr(var); A=detr(a);
                if (isa(var{1},'numeric') && isa(A{1},'numeric'))
                    out{jj,ii}=(-1)^(ii+jj)*var{1}/A{1};
                else
                    ss=Opr({(-1)^(ii+jj)},Opr(var,A,'/'),'*');
                    out{jj,ii}=Modify(ss{1});
                end
            end
        end
    case 'n'
        out=[];
        if (n1~=1)
            error('Not defined');
        end
        count_ch=[];
        num=0; ch=[];
        for ii=1:m1
            if (isa(a{ii},'char'))
                count_ch=[count_ch,ii];
                if (~isempty(ch))
                    ch=sprintf('%s+%s^2',ch,Parantes2(a{ii},'n'));
                else
                    ch=sprintf('%s^2',Parantes2(a{ii},'n'));
                end
            elseif (isa(a{ii},'numeric'))
                %                 count_num=[count_num,ii];
                num=num+a{ii}^2;
            end
            if (isempty(count_ch))
                out=sqrt(num);
            elseif (num==0)
                out=sprintf('sqrt(%s)',ch);
            else
                out=sprintf('sqrt(%d+%s)',num,ch);
            end
        end
        %         out=num2cell(out);
    case '/'
        if (size(b)~=[1 1])
            error('Not defined');
        end
        for ii=1:m1
            for jj=1:n1
                if (isa(a{ii,jj},'numeric') && isa(b{1},'numeric'))
                    out{ii,jj}=a{ii,jj}/b{1};
                else
                    if v2s(b{1}==1)
                        %                         out{ii,jj}=sprintf('%s',Parantes(num2str(a{ii,jj})));
                        out{ii,jj}=sprintf('%s',num2str(a{ii,jj}));
                    elseif v2s(b{1}==-1)
                        out{ii,jj}=Modify(sprintf('-%s',Parantes(num2str(a{ii,jj}))));
                    elseif v2s(b{1}==0)
                        error('Not defined num/0');
                    elseif v2s(a{ii,jj}==1)
                        out{ii,jj}=sprintf('1/%s',Parantes2(num2str(b{1}),'/'));
                    elseif v2s(a{ii,jj}==-1)
                        out{ii,jj}=sprintf('1/%s',Modify(sprintf('-%s',Parantes2(num2str(b{1}),'/'))));
                    elseif v2s(a{ii,jj}==0)
                        out{ii,jj}=0;
                    else
                        out{ii,jj}=sprintf('%s/%s',Parantes(num2str(a{ii,jj})),Parantes2(num2str(b{1}),'/'));
                    end
                end
            end
        end
    case 'c'
        if (n1~=1 || n2~=1 || m1~=m2)
            error('Not defined');
        end
        var=Opr(Opr({a{2}},{b{3}},'*'),Opr({a{3}},{b{2}},'*'),'-');
        out{1,1}=var{1};
        var=Opr(Opr({a{3}},{b{1}},'*'),Opr({a{1}},{b{3}},'*'),'-');
        out{2,1}=var{1};
        var=Opr(Opr({a{1}},{b{2}},'*'),Opr({a{2}},{b{1}},'*'),'-');
        out{3,1}=var{1};
    case '^'
        if (isa(a{1},'char'))
            out{1}=sprintf('%s^%d',Parantes2(a{1},'^'),b);
        else
            out{1}=a{1}^b;
        end
end