function [Covariances]=Covariance_init(Clusters,Means,i,u,Covariances)
% Covariance intialization considering if the final matrix will be singular

 % Standard method
% I method (Andrew/EM-init)
%     for j=1:length(Clusters.(u))            
%         P=(Clusters.(u)(j,:)-Means(i,:))'*(Clusters.(u)(j,:)-Means(i,:));
%         Covariances(:,:,i)=(Covariances(:,:,i)+P);
%     end
%     Covariances(:,:,i)= (1/(length(Clusters.(u))-1))*Covariances(:,:,i);

    Covariances(:,:,i) = cov(Clusters.(u)); 

% Cov is close to singular-> Cov. = spherical cov
%      [~,r] = chol(Covariances(:,:,i));
%       r == 0 && rank(Covariances(:,:,i)) == size(Covariances(:,:,i),1)
%      else

% I Method 

%         DIM=length(Means);
%         cova = eye(DIM,DIM);
%         SS = 0;
% 
%                  for k=1:DIM
%                      for j=1:length(Clusters.(u))
%                      S=(Clusters.(u)(j,k)-Means(i,k)).^2;
%                      SS=SS+S;
%                      end
%                      cova(k,k)=(1/((length(Clusters.(u)-1))))*SS;
%                      SS=0;
%                  end
%                  Covariances(:,:,i)=cova;
         
% II METHOD (paper EM-init)
%           for j=1:length(Clusters.(u))
%                  S=sum((Clusters.(u)(j,:)-Means(i,:)).^2)*(eye(size(Clusters.(u),2),size(Clusters.(u),2)));
%                  Covariances(:,:,i)=(Covariances(:,:,i)+S);
%           end
%           Covariances(:,:,i)=(1/(length(Clusters.(u))*length(Covariances(:,:,i))))*Covariances(:,:,i);
%      end
%      
% Even the Sphe.covariance gives a Singular matrix -> Cov = identity matrix
% I METHOD
%      [~,r] = chol(Covariances(:,:,1));
%      if r == 0 && rank(Covariances(:,:,1)) == size(Covariances(:,:,1),1)
%      else
%         Covariances(:,:,i)=eye(size(Clusters.(u),2),size(Clusters.(u),2));
%      end
     
    
end