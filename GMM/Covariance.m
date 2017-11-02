function [Covariances_new]=Covariance(Training_set,Means_new,Covariances_new,j,Weights)
% Covariance intialization considering if the final matrix will be singular

     % Standard method
%     I METHOD (Adndrew)
%      for i = 1:length(Training_set)
%          P =(Training_set(i,:)-Means_new(:,j)')'*(Training_set(i,:)-Means_new(:,j)');
%          Covariances_new(:,:,j) = Covariances_new(:,:,j)+P;
%      end
%      Covariances_new(:,:,j)= (1/((size(Training_set,1))-1))*Covariances_new(:,:,j);

%     II METHOD (video Victor Lavrenko)
     for i = 1:length(Training_set)
         P = Weights(i,j)*(Training_set(i,:)-Means_new(:,j)')'*(Training_set(i,:)-Means_new(:,j)');
         Covariances_new(:,:,j) = Covariances_new(:,:,j)+P;
     end
%      Covariances_new(:,:,j)= (1/(sum(Weights(:,j))))*Covariances_new(:,:,j);
%      
 % Cov is close to singular-> Cov. = spherical cov
%      [~,r] = chol(Covariances_new(:,:,j));
%      if r == 0 && rank(Covariances_new(:,:,j)) == size(Covariances_new(:,:,j),1)
%      else
% I METHOD
%          for i=1:length(Training_set)
%              S=sum(Weights(i,j)*(Training_set(i,:)-Means_new(:,j)').^2)*(eye(size(Training_set,2),size(Training_set,2)));
%              Covariances_new(:,:,j)=(Covariances_new(:,:,j)+S);
%          end
%              Covariances_new(:,:,j)=(1/(length(Training_set)*length(Covariances_new(:,:,j))))*Covariances_new(:,:,j);
% II METHOD

%     DIM=length(Means_new);
%     cova = eye(DIM,DIM);
%     SS = 0;
% 
%              for k=1:DIM
%                  for j=1:length(Training_set)
%                  S = Weights(j,q)*(Training_set(j,k)-Means_new(k,1)).^2;
%                  SS=SS+S;
%                  end
%                  cova(k,k)=(1/((length(Training_set-1))))*SS;
%                  SS=0;
%              end
% 
%              Covariances_new(:,:,q)=cova;
         
%      end

% Even the Sphe.covariance gives a Singular matrix -> Cov = identity matrix

%      [~,r] = chol(Covariances_new(:,:,j));
%      if r == 0 && rank(Covariances_new(:,:,j)) == size(Covariances_new(:,:,j),1)
%      else
%          Covariances_new(:,:,j)=eye(size(Training_set,2),size(Training_set,2));
%      end
%      

  
end