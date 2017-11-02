function    [V_new,Prob_new]=Viterbi_Composite2(HMM,HMMs,Data,i,Gesture,K,V_old,Prob_old)

if i==1    
    for j=1:Gesture
        if j ~= 7
            clear Pb;
            Training_set = Data;
            Covariances = HMMs(1,j).Covariances;
            Means = HMMs(1,j).Mean;
            for k=1:K
                Pb(:,k)= mvnpdf(Training_set,Means(k,:),Covariances(:,:,k));     
            end          
            v(j,1)=HMM.Beginning_P(1,j)*sum(HMMs(1,j).Beginning_P'.*HMMs(1,j).Prior.*Pb);
        else
            v(j,1)=0;
        end
    end
    [Prob_new,V_new]=max(v);
else
     for j=1:Gesture
        if j ~= 7 
            clear Pb;
            Training_set = Data;
            Covariances = HMMs(1,j).Covariances;
            Means = HMMs(1,j).Mean;
            for k=1:K
                Pb(:,k)= mvnpdf(Training_set,Means(k,:),Covariances(:,:,k));     
            end
            if j == V_old
                tras = (HMMs(1,j).Trans_P.*HMMs(1,j).Prior.*Pb);
                SUM(j,1)=sum(sum((tras)));
            else
                SUM(j,1)=sum(HMMs(1,j).Beginning_P'.*HMMs(1,j).Prior.*Pb);
            end
        else
            SUM(j,1)=0;
        end
     end  
     
    v=Prob_old.*HMM.Trans_P(V_old,:).*SUM';
    [Prob_new,V_new]=max(v);

end

end

