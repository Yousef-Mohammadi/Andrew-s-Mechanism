% function [ Q,Qdot ] = Andrew_system( x,step,t )
% clear all;
clc;
disp('individual');
close all;
global n_d
n_q=1; n_d=7; n_j=7; g=9.81; 
% AG=0.04; GE=0.02; AH=0.04; HE=0.02; BE=0.035; EF=0.028; OF=0.007; C_0=4530;
% OA_x=-.06934; OA_y=-.00227; OB_x=-.03635; OB_y=.03273; OC_x=.014; OC_y=.072;
sc=0.018; sd=0.02; xc=0.014; yc=0.072;
param={'AG','GE','OA_x','OA_y','OB_x','OB_y','OC_x','OC_y','AH','HE','BE','EF','OF'};
% angle={'th1','th2','th3','th4','th5','th6','th7'};
param_eval={0.04,0.02,-0.06934,-0.00227,-0.03635,0.03273,0.014,0.072,0.04,0.02,0.035,0.028,0.007};
PARAM=[param;param_eval];
G={'OA_x+AH*cos(th2)+HE*cos(th2+th6)+EF*cos(th2+th6+th7)-OF*cos(th1)';
    'OA_y+AH*sin(th2)+HE*sin(th2+th6)+EF*sin(th2+th6+th7)-OF*sin(th1)';
    'AG*cos(th3)+GE*cos(th3+th5)-AH*cos(th2)-HE*cos(th2+th6)';
    'AG*sin(th3)+GE*sin(th3+th5)-AH*sin(th2)-HE*sin(th2+th6)';
    'OA_x+AH*cos(th2)+HE*cos(th2+th6)-BE*cos(th4)-OB_x';
    'OA_y+AH*sin(th2)+HE*sin(th2+th6)-BE*sin(th4)-OB_y';
    'q-th1'};
G=sym(G);
for ii=1:length(G)
    eval(sprintf('G(ii)=subs(G(ii),''%s'',%f); ',PARAM{:}));
end
eps1=eye(n_d)*1e-6; eps2=eye(n_q)*1e-6;  eps=1e-6;
Q=[]; Qdot=[]; Q(1,1)=0.0617; Qdot(1,1)=0; q_eval=0.0617; qdot_eval=0; Ang=[];
eval(sprintf('%s=%f;',PARAM{:}));
% options=optimoptions('fsolve','TolFun',1e-13,'TolX',1e-10,'MaxFunEvals',500*n_d,'MaxIter',1000);
options=optimoptions('fmincon','Algorithm','sqp','TolCon',1e-4,'MaxFunEvals',100*n_d,'MaxIter',1000);
cc=0; color=colormap;
% step=.05/100;
%% gradiant projection:
for tt=t
    fprintf('time=%d\n',tt);
    for ll=1:4
        flag=-2;
%         eval(sprintf('G_s(%d)=subs(G_s(%d),''q'',q_eval);',[1:length(G);1:length(G)]));
%         Sol=solve(G_s,sym(angle));
        while flag~=1
%             if tt==0
%                 init=randn(1,7);
%             else
%                 init=Ang(end,:);
%             end
            try
                [Sol,~,flag]=fmincon(@(x)0,randn(7,1),[],[],[],[],[],[],@(x)NLC(x,q_eval),options);
%                 [Sol,~,flag]=fsolve(@(x)NLC(x,q_eval),rand(7,1));
            catch
                flag=-2;
            end
        end
        [~,fval]=NLC(Sol,q_eval);
        ss=sprintf('%d,',fval);
        fprintf('fval=[%s] && exitflag=%i\n',ss,flag);
        if flag~=1
            error('flag~=1');
        end
%         ss=[num2cell(1:n_j);angle];
        ss=[num2cell(1:n_j);num2cell(1:n_j)];
%         eval(sprintf('angle_eval(%d)=double(Sol.%s);',ss{:}));
        eval(sprintf('angle_eval(%d)=double(Sol(%d));',ss{:}));
        if ll==1
            Ang=[Ang;angle_eval];
        end
        eval(sprintf('th(%d)=angle_eval(%d);',[1:n_d;1:n_d]));
%         ghat=zeros(n_d,n_d); dghattoq=zeros(n_d*n_q,n_q);
        for ii=1:n_d
            for jj=1:n_d
                var1=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval]);
                var2=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                ghat(ii,jj)=(var1-var2)/eps;
%                 for kk=1:n_q
%                     var1=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval+eps2(kk,:)]);
%                     var2=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval]);
%                     var3=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(kk,:)]);
%                     var4=Runge_Kutta(G(ii),[angle_eval,q_eval]);
%                     dghattoq((ii-1)*n_q+kk,jj)=(var1-var2-var3+var4)/eps^2;
%                 end
            end
        end
        dghattoq=zeros(7,7);
%         dGtoq=zeros(n_d*n_q,1); ddGtoq=zeros(n_d*n_q^2,1);
        for ii=1:n_d
            for jj=1:n_q
                var1=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(jj,:)]);
                var2=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                dGtoq((ii-1)*n_q+jj,1)=(var1-var2)/eps;
                for kk=1:n_q
                    var1=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(jj,:)+eps2(kk,:)]);
                    var2=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(jj,:)]);
                    var3=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(kk,:)]);
                    var4=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                    ddGtoq((ii-1)*n_q^2+(jj-1)*n_q+kk,1)=(var1-var2-var3+var4)/eps^2;
                end
            end
        end
        gradf=Prediction(ghat,dGtoq,n_q);
        var=zeros(n_d*n_q^2,1);
        for ii=1:n_d  
            var=var+Diadic(eye(n_d),FU(n_q,n_q))*Diadic(dghattoq(:,ii),gradf(:,ii));
        end
        hessian=Prediction(ghat,ddGtoq+var,n_q^2);
        %% Transformations: 
        Len={3,[3,5],2,[2,6],4,[2,6,7],1};
        [A0toj,Ajtoi,Ahatj]=Transform(th,PARAM);
        %% DBtoq:
        for ii=1:n_d
            for jj=length(Len{ii})
                dbtoq=zeros(4*n_q,4);
                for kk=1:jj
                    W(:,:,kk,jj,ii)=A0toj(:,:,kk,ii)*Ahatj(:,:,kk)*Ajtoi(:,:,kk,jj,ii);
                    dbtoq=dbtoq+Diadic(W(:,:,kk,jj,ii)',gradf(:,Len{ii}(kk)));
                end
            dBtoq(:,:,jj,ii)=dbtoq;
            end
        end
        %% DDBtoq:
        for ii=1:n_d
            for kk=length(Len{ii})
%                 var=zeros(4*n_q,4);
                ddbtoq=zeros(4*n_q^2,4);
                for jj=1:kk
                    var=zeros(4*n_q,4);
                    ddbtoq=ddbtoq+Diadic((W(:,:,jj,kk,ii))',hessian(:,Len{ii}(jj)));
                    for ss=1:kk  
                        if (ss<=jj)
                            What=W(:,:,ss,jj,ii)*Ahatj(:,:,jj)*Ajtoi(:,:,jj,kk,ii);
                            var=var+Diadic(What',gradf(:,Len{ii}(jj)));
                        elseif (ss>jj)
                            What=W(:,:,jj,ss,ii)*Ahatj(:,:,ss)*Ajtoi(:,:,ss,kk,ii);
                            var=var+Diadic(What',gradf(:,Len{ii}(ss)));
                        end
%                         ddbtoq=ddbtoq+Diadic(What',Diadic(gradf(:,Len{ii}(jj)),gradf(:,Len{ii}(ss))));
                    end
                    ddbtoq=ddbtoq+Diadic(eye(4),FU(n_q,n_q))*Diadic(var,gradf(:,Len{ii}(jj)));
                end
                ddBtoq(:,:,kk,ii)=ddbtoq;
            end
        end
        m_c=[0.0705,0.00706,0.05498,0.00706,0.02373,0.00365,0.04325];
        I(:,:,1)=[0,0,0;0,0,0;0,0,1.169e-5]; I(:,:,2)=[0,0,0;0,0,0;0,0,5.667e-7];
        I(:,:,3)=[0,0,0;0,0,0;0,0,1.912e-5]; I(:,:,4)=[0,0,0;0,0,0;0,0,5.667e-7];
        I(:,:,5)=[0,0,0;0,0,0;0,0,5.255e-6]; I(:,:,6)=[0,0,0;0,0,0;0,0,4.41e-7];
        I(:,:,7)=[0,0,0;0,0,0;0,0,2.194e-6]; 
        rc_c=[[0.02308;0.00916;0],[0.0058;0;0],[0.01228;-0.00449;0],[0.0058;0;0],[0.01874;0.01043;0],[0.0165;0;0],[0.00092;0;0]];
        for ii=1:n_j
            val=rc_c(:,ii)'*rc_c(:,ii);
            Ip(:,:,ii)=I(:,:,ii)+m_c(ii)*(val*eye(3)-rc_c(:,ii)*rc_c(:,ii)');
            Ipp(:,:,ii)=-(Ip(:,:,ii)-0.5*trace(Ip(:,:,ii))*eye(3));
        end
        C1=Diadic(Vec((eye(4))),(eye(n_q)));
        C2=Diadic((eye(4)),FU(n_q,4))*C1;
        C3=Diadic(Diadic(eye(4),qdot_eval),FU(n_q,4))*C1;
        for ii=1:n_j
            J(:,:,ii)=[Ipp(:,:,ii),m_c(ii)*rc_c(:,ii);m_c(ii)*rc_c(:,ii)',m_c(ii)];
            for jj=length(Len{ii})
                M(:,:,jj,ii)=C1'*Diadic(J(:,:,ii),FU(n_q,4)')*Diadic(dBtoq(:,:,jj,ii)*dBtoq(:,:,jj,ii)',eye(4))*C2;
                N(:,:,jj,ii)=C1'*Diadic(J(:,:,ii),FU(n_q,4)')*Diadic(dBtoq(:,:,jj,ii)*ddBtoq(:,:,jj,ii)',eye(4))*C3;
                GG(:,jj,ii)=m_c(ii)*Diadic([rc_c(:,ii)',1],eye(n_q))*dBtoq(:,:,jj,ii)*[0;-g;0;0];
            end
        end
        N_total=zeros(n_q,n_q);
        M_total=zeros(n_q,n_q);
        GG_total=zeros(n_q,1);
        for ii=1:n_j
            finish=length(Len{ii});
%             for jj=1:finish
%                 M_total=M_total+M(:,:,jj,ii);
                M_total=M_total+M(:,:,finish,ii);
%                 N_total=N_total+N(:,:,jj,ii);
                N_total=N_total+N(:,:,finish,ii);
%                 GG_total=GG_total+GG(:,jj,ii);
                GG_total=GG_total+GG(:,finish,ii);
%             end
        end
        spring_org=[xc;yc;0];
        spring_ins=[sc;sd;0];
%         fun=1;
%         eval(sprintf('fun=[fun,sin(%d*pi*tt/t(end))];',1:10));
%         eval(sprintf('fun=[fun,cos(%d*pi*tt/t(end))];',1:10));
%         eval(sprintf('fun=[fun,heaviside(tt-%f-1e-5)];',t));
        l0=0.07785; c0=4530; M0=0.033;
        beta1=gradf(:,4)*(A0toj(1:3,1:3,1,5)'*[0;0;1])';
        beta2=gradf(:,1)*(A0toj(1:3,1:3,1,7)'*[0;0;1])';
        var=A0toj(:,:,1,5)*[spring_ins;1];
%         OB=[OB_x;OB_y;0];
        rp_F=var(1:3);
        rp_L=spring_org;
        len=norm(rp_F-rp_L);
        LA=(rp_L-rp_F)/len;
        Ehat1=cross(rp_F,LA);
        Ehat2=[0;0;1];
        R_spring=beta1*Ehat1;
        F=-c0*(len-l0);
        T_ext1=R_spring*F;
        T_ext2=beta2*Ehat2*M0;
        fprintf('T_ext1=%d & T_ext2=%d & F=%d and len-l0=%d\n',T_ext1,T_ext2,F,len-l0);
        if abs(T_ext1)>abs(T_ext2)
%             pause;
        end
        T_ext=T_ext1+T_ext2;
        %% Forward Dynamic
%         figure(1);
        if ll==1
            cc=cc+1;
            fprintf('cc=%d\n',cc);
            PlotStructure(A0toj,OF,cc);
            title(sprintf('Time=%f',tt));
            xlabel('X');
            ylabel('Y');
            pause(0.01);
            cla
            K1_q=qdot_eval*step; K1_qd=(M_total\(T_ext+GG_total-N_total*qdot_eval))*step;
%             q_eval=Q(fix(tt/step)+1,:)'+K1_q/2; qdot_eval=Qdot(fix(tt/step)+1,:)'+K1_qd/2;
            q_eval=Q(end,:)'+K1_q/2; qdot_eval=Qdot(end,:)'+K1_qd/2;
        elseif ll==2
            K2_q=qdot_eval*step; K2_qd=(M_total\(T_ext+GG_total-N_total*qdot_eval))*step;
            q_eval=Q(end,:)'+K2_q/2; qdot_eval=Qdot(end,:)'+K2_qd/2;
        elseif ll==3
            K3_q=qdot_eval*step; K3_qd=(M_total\(T_ext+GG_total-N_total*qdot_eval))*step;
            q_eval=Q(end,:)'+K3_q; qdot_eval=Qdot(end,:)'+K3_qd;
        else
            K4_q=qdot_eval*step; K4_qd=(M_total\(T_ext+GG_total-N_total*qdot_eval))*step;
        end
    end
    q_eval=Q(end,:)'+(K1_q+2*K2_q+2*K3_q+K4_q)/6;
    qdot_eval=Qdot(end,:)'+(K1_qd+2*K2_qd+2*K3_qd+K4_qd)/6;
    if q_eval>pi
        q_eval=q_eval-2*pi;
    end
    Q(end+1,:)=q_eval; Qdot(end+1,:)=qdot_eval;
end