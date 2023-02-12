for ii=1:n_d
            for jj=1:n_d
                var1=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval]);
                var2=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                ghat(ii,jj)=(var1-var2)/eps;
                for kk=1:n_q
                    var1=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval+eps2(kk,:)]);
                    var2=Runge_Kutta(G(ii),[angle_eval+eps1(jj,:),q_eval]);
                    var3=Runge_Kutta(G(ii),[angle_eval,q_eval+eps2(kk,:)]);
                    var4=Runge_Kutta(G(ii),[angle_eval,q_eval]);
                    dghattoq((ii-1)*n_q+kk,jj)=(var1-var2-var3+var4)/eps^2;
                end
            end
        end
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