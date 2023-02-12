function [ C ] = Prediction( A,B,n_Q )
global n_d
gama(:,:,n_d+1)=eye(n_d);
for ii=n_d:-1:2
    sigma=double(eye(n_d)-A(:,ii)*A(:,ii)'*gama(:,:,ii+1)/norm(gama(:,:,ii+1)*A(:,ii))^2);
    gama(:,:,ii)=double(gama(:,:,ii+1)*sigma);
end
for ii=1:n_d
    if (ii==1)
        C(:,ii)=double(-Diadic(A(:,ii)'*gama(:,:,ii+1),eye(n_Q))*B/norm(gama(:,:,ii+1)*A(:,ii))^2);
        eta(:,ii)=double(Diadic(A(:,ii),C(:,ii)));
    else
        C(:,ii)=double(-Diadic(A(:,ii)'*gama(:,:,ii+1),eye(n_Q))*(B+eta(:,ii-1))/norm(gama(:,:,ii+1)*A(:,ii))^2);
        eta(:,ii)=double(eta(:,ii-1)+Diadic(A(:,ii),C(:,ii)));
    end
end