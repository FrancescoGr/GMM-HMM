function [obj]= GMM_modeling_built_new(Training_set)

%%  Initialization

% 1) Number of Gaussians in the model
K = 3;
it = 1500;
Reg = 0;

        SharedCovariance = {false};
        options = statset('MaxIter',1000,'Display','final'); % Increase number of EM iterations
        obj = fitgmdist(Training_set,K,'CovarianceType','full',...
            'SharedCovariance',false,'Options',options,'Regularize',Reg,'Replicates', it);
        

end