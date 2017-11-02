function [DatawithStates,Emission_P_Sequence]=Labeling(GMMs,Gesture,States)
for i=1:lenght(Gesture)
    %% Initializations
    Covariances=GMMs(1,i).Covariances;
    Mean=GMMs(1,i).Mean;
    Priors=GMMs(1,i).Priors;
    K=States(1,i);
    
    Data_set=Gesture{1,i};
    LD=length(Data_set);
    
    %% Labeling  
    
    Pb=[];
            % Probability of Belonging in the K state
            for j = 1:K
                for i=1:LD
                    Pb(i,j) = (1/((sqrt((2*pi)^(n)*det(Covariances(:,:,j))))))*exp(-0.5*(Data_set(i,:)-Means(j,:))*pinv(Covariances(:,:,j))*(Data_set(i,:)-Means(j,:))');
                end
            end

            Pp = [];
            index_class = [];
            for i = 1:LD
                Sum = 0;
                for k = 1:K
                    Sum = Sum+(Pb(i,k)*Priors(k,1));
                end
                % Posterior probability
                for j = 1:K
                    Pp(i,j) = (Pb(i,j)*Priors(j,1))/Sum;
                end 


                % State prediction for the i-sample
                [~,index_class(i,1)]=max(Pp(i,:));
            end
            
%             Weights=[];
%             for j = 1:K
%                 for i = 1:LD
%                         Weights(i,j) = Pp(i,j)/sum(Pp(:,j));
%                 end   
%             end

         DatawithStates{1,i} = [Data_set, index_class ];
         Emission_P_Sequence{1,i}(:,:)=Pb*Priors;
         
end
end