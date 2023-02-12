function  [C,Ceq]  = NLC( th,q )
global n_d
    param={'AG','GE','OA_x','OA_y','OB_x','OB_y','OC_x','OC_y','AH','HE','BE','EF','OF'};
    param_eval={0.04,0.02,-0.06934,-0.00227,-0.03635,0.03273,0.014,0.072,0.04,0.02,0.035,0.028,0.007};
    PARAM=[param;param_eval];
    eval(sprintf('%s=%d;',PARAM{:}));
    G={'OA_x+AH*cos(th(2))+HE*cos(th(2)+th(6))+EF*cos(th(2)+th(6)+th(7))-OF*cos(th(1))';
        'OA_y+AH*sin(th(2))+HE*sin(th(2)+th(6))+EF*sin(th(2)+th(6)+th(7))-OF*sin(th(1))';
        'AG*cos(th(3))+GE*cos(th(3)+th(5))-AH*cos(th(2))-HE*cos(th(2)+th(6))';
        'AG*sin(th(3))+GE*sin(th(3)+th(5))-AH*sin(th(2))-HE*sin(th(2)+th(6))';
        'OA_x+AH*cos(th(2))+HE*cos(th(2)+th(6))-BE*cos(th(4))-OB_x';
        'OA_y+AH*sin(th(2))+HE*sin(th(2)+th(6))-BE*sin(th(4))-OB_y';
        'q-th(1)'};
    for ii=1:length(G)
        F(ii)=eval(G{ii});
    end
    Ceq=F; 
    C=1e-1-abs(th(2)-th(3)); C=[C;abs(th(2)-th(3))-2*pi+1e-1];
%     C=[C;th(2);-pi/2-th(2)];
%     C=[C;-th(3);-pi/2+th(3)];
%     C=[];
%     eval(sprintf('C=[C;th(%d)-pi-3.5e-1];',1));
%     eval(sprintf('C=[C;-th(%d)-pi-3.5e-1];',1));
%     eval(sprintf('C=[C;th(%d)-2*pi-1e-1];',1:n_d));
%     eval(sprintf('C=[C;-th(%d)-2*pi-2e-1];',2:n_d));
end
