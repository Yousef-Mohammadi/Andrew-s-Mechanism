function [ A0toj,Ajtoi,Ahatj ] = Transform( th,PARAM )
    global n_d
    eval(sprintf('%s=%f;',PARAM{:}));
    eval(sprintf('th%d=th(%d);',[1:n_d;1:n_d]));
    A0toj(:,:,1,1)=[cos(th3),-sin(th3),0,OA_x;sin(th3),cos(th3),0,OA_y;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,1)=eye(4);
    Ahatj(:,:,1)=[0,-sign(th3),0,0;sign(th3),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,2)=A0toj(:,:,1,1);
    Ajtoi(:,:,1,2,2)=[cos(th5),-sin(th5),0,AG;sin(th5),cos(th5),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,2,2)=A0toj(:,:,1,2)*Ajtoi(:,:,1,2,2);
    Ajtoi(:,:,2,2,2)=eye(4);
    Ajtoi(:,:,1,1,2)=eye(4);
    Ahatj(:,:,2)=[0,-sign(th5),0,0;sign(th5),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,3)=[cos(th2),-sin(th2),0,OA_x;sin(th2),cos(th2),0,OA_y;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,3)=eye(4);
    Ahatj(:,:,3)=[0,-sign(th2),0,0;sign(th2),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,4)=A0toj(:,:,1,3);
    Ajtoi(:,:,1,2,4)=[cos(th6),-sin(th6),0,AH;sin(th6),cos(th6),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,2,4)=A0toj(:,:,1,4)*Ajtoi(:,:,1,2,4);
    Ajtoi(:,:,1,1,4)=eye(4);
    Ajtoi(:,:,2,2,4)=eye(4);
    Ahatj(:,:,4)=[0,-sign(th6),0,0;sign(th6),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,5)=[cos(th4),-sin(th4),0,OB_x;sin(th4),cos(th4),...
        0,OB_y;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,5)=eye(4);
    Ahatj(:,:,5)=[0,-sign(th4),0,0;sign(th4),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,6)=A0toj(:,:,1,3); A0toj(:,:,2,6)=A0toj(:,:,2,4); 
    Ajtoi(:,:,2,3,6)=[cos(th7),-sin(th7),0,HE;sin(th7),cos(th7),...
        0,0;0,0,1,0;0,0,0,1];
    A0toj(:,:,3,6)=A0toj(:,:,2,6)*Ajtoi(:,:,2,3,6);
    Ajtoi(:,:,1,1,6)=eye(4);
    Ajtoi(:,:,2,2,6)=eye(4); 
    Ajtoi(:,:,1,2,6)=Ajtoi(:,:,1,2,4);
    Ajtoi(:,:,3,3,6)=eye(4);
    Ajtoi(:,:,1,3,6)=Ajtoi(:,:,1,2,4)*Ajtoi(:,:,2,3,6);
    Ahatj(:,:,6)=[0,-sign(th7),0,0;sign(th7),0,0,0;0,0,0,0;0,0,0,0];
    A0toj(:,:,1,7)=[cos(th1),-sin(th1),0,0;sin(th1),cos(th1),0,0;0,0,1,0;0,0,0,1];
    Ajtoi(:,:,1,1,7)=eye(4);
    Ahatj(:,:,7)=[0,-sign(th1),0,0;sign(th1),0,0,0;0,0,0,0;0,0,0,0];
%     Ahatj=[0,-1,0,0;1,0,0,0;0,0,0,0;0,0,0,0];
end

