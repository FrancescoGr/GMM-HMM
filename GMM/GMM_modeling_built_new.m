function [obj]= GMM_modeling_built_new(Training_set,choice)

%%  Initialization

% 1) Number of Gaussians in the model
if choice==1
    it = 1700;
else
    if size(Training_set,1)<10000
        it = 15000;
    else
        it = 10000;
    end
end

K = 3;
Reg = 0;

        SharedCovariance = {false};
        options = statset('MaxIter',1000,'Display','final'); % Increase number of EM iterations
        obj = fitgmdist(Training_set,K,'CovarianceType','full',...
            'SharedCovariance',false,'Options',options,'Regularize',Reg,'Replicates'
        );
        

end