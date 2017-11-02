function [Clusters]= Partitions(Training_set,K)
% Defines the first partition of the dataset in K clusters

% Initializations
for k=1:K
    u = ['Clu',int2str(k)];
    Clusters.(u)=[];
end

Data_length=length(Training_set);

% Generation of K rand points defined as centroids of the clusters
Check =1;
while Check ~= 0
    Random=randi([1,Data_length],K,1);
    Check=length(find(Random==0));
end

% Computation of the distances from centroids
for i=1:Data_length
    for j=1:K
        dist(j,1)=sum(Training_set(i,:)-Training_set(Random(j,1),:));
    end
    
        
        n_clu = find(abs(dist)== min(abs(dist)));
        n_clu = n_clu(1,1);
        u = ['Clu',int2str(n_clu)];
    
    % Cluster division
    Clusters.(u)=[Clusters.(u); Training_set(i,:)];
end
end