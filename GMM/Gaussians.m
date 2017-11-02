function [Covariances,Means,Priors]=Gaussians(K,Training_set,Clusters)
% Initialization of the K Gaussians parameters

Covariances=zeros(size(Training_set,2),size(Training_set,2),K);

for i=1:K
    u=['Clu',int2str(i)];
    
    % Means definition
    Means(i,:)= mean(Clusters.(u));
    
    % Priors definition
    Priors(i,1)=length(Clusters.(u))/length(Training_set);
    
    % Covariances
    [Covariances]=Covariance_init(Clusters,Means,i,u,Covariances);

%     DIM=length(Means);
%     cova = eye(DIM,DIM);
%     SS = 0;
% 
%          for k=1:DIM
%              for j=1:length(Training_set)
%              S = (Training_set(j,k)-Means(i,k)).^2;
%              SS=SS+S;
%              end
%              cova(k,k)=(1/(DIM*(length(Training_set-1))))*SS;
%              SS=0;
%          end
%          
%          Covariances(:,:,i)=cova;

    
end
end