function    [V_new,Prob_new,MultiFactor_new]=Viterbi_Composite(HMM,HMMs,Data,MultiFactor_old,i,Gesture,K,V_old,Prob_old,Sum_Pb)
if i==1    
    for j=1:Gesture
        if j ~= 7
            clear Pb;
            for k=1:K
                Training_set = Data;
                Covariances = HMMs(1,j).Covariances;
                Means = HMMs(1,j).Mean;
                n = size(Means,2);
                Pb(1,k) = (1/((sqrt((2*pi)^(n)*det(Covariances(:,:,k))))))*exp(-0.5*(Training_set-Means(k,:))*pinv(Covariances(:,:,k))*(Training_set-Means(k,:))');
            end
            Pb=Pb./Sum_Pb;           
            v(j,1)=HMM.Beginning_P(1,j)*sum(HMMs(1,j).Beginning_P'.*HMMs(1,j).Prior.*Pb);
            MultiFactor(j,:)=HMMs(1,j).Beginning_P'.*HMMs(1,j).Prior.*Pb;
        else
            v(j,1)=0;
        end
    end
    [Prob_new,V_new]=max(v);
%     MultiFactor_new = HMMs(1,V_new).Beginning_P;
    MultiFactor_new = MultiFactor(V_new,:);
else
     for j=1:Gesture
        if j ~= 7 
            clear Pb;
            for k=1:K
                Training_set = Data;
                Covariances = HMMs(1,j).Covariances;
                Means = HMMs(1,j).Mean;
                n = size(Means,2);
                Pb(1,k) = (1/((sqrt((2*pi)^(n)*det(Covariances(:,:,k))))))*exp(-0.5*(Training_set-Means(k,:))*pinv(Covariances(:,:,k))*(Training_set-Means(k,:))');
            end
            
            Pb=Pb./Sum_Pb;

            if j == V_old
                tras = (HMMs(1,j).Trans_P.*HMMs(1,j).Prior.*Pb);
                tras1= sum(tras.*MultiFactor_old);
                SUM(j,1)=sum(sum((tras.*MultiFactor_old)'));
            else
                SUM(j,1)=sum(HMMs(1,j).Beginning_P'.*HMMs(1,j).Prior.*Pb);
                MultiFactor(j,:)=HMMs(1,j).Beginning_P'.*HMMs(1,j).Prior.*Pb;
            end
        else
            SUM(j,1)=0;
        end
     end  
     
    v=Prob_old.*HMM.Trans_P(V_old,:).*SUM';
    [Prob_new,V_new]=max(v);
    
    if V_new==V_old
        MultiFactor_new = tras1';
    else
        MultiFactor_new = MultiFactor(V_new,:); % HMMs(1,V_new).Beginning_P;
    end
end

end