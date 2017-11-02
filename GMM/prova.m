clear 
rng(0,'twister')  % For reproducibility
mu1 = [1 2];
sigma1 = [3 .2; .2 2];
mu2 = [-1 -2];
sigma2 = [2 0; 0 1];
X = [mvnrnd(mu1,sigma1,200); mvnrnd(mu2,sigma2,100)];
gm = fitgmdist(X,2);
threshold = [0.4 0.6];
P = posterior(gm,X);
gmSharedDiag = fitgmdist(X,2,'CovType','Diagonal','SharedCovariance',true');
[idxSharedDiag,~,PSharedDiag] = cluster(gmSharedDiag,X);




%% to have obj
for i=1:15
SharedCovariance = {false};
options = statset('MaxIter',1000,'Display','final'); % Increase number of EM iterations
obj{i,1} = fitgmdist(X,2,'CovarianceType','full',...
            'SharedCovariance',false,'Options',options,'Regularize',0.0000000001,'Replicates', 500);
  
end

        %% to have correctly log-likehood
covNames = { 'diagonal','full'};
CovType = find(strncmpi(GMMs.CovType,covNames,length(GMMs.CovType)));

    mu=GMMs.mu;
    log_prior = log(GMMs.PComponents);
    [n,d]=size(X);
    k=size(GMMs.mu,1);
    Sigma=GMMs.Sigma;
    
    logDetSigma = -Inf;
    
    for j = 1:k
                % compute the log determinant of covariance
                [L,f] = chol(Sigma(:,:,j) );
                diagL = diag(L);
                if (f ~= 0) || any(abs(diagL) < eps(max(abs(diagL)))*size(L,1))
                     error(message('stats:gmdistribution:wdensity:IllCondCov'));
                end
                logDetSigma = 2*sum(log(diagL));
                log_lh(:,j) = sum(((X - mu(j,:))/L).^2, 2); 
                log_lh(:,j) = -0.5*(log_lh(:,j) + logDetSigma);
    end
                log_lh = log_lh + log_prior - d*log(2*pi)/2;
                
                
maxll = max(log_lh,[],2);
%minus maxll to avoid underflow
% post = exp(log_lh-maxll);
Pb_correct=exp(log_lh);

        %% to have correctly Pb
        
for j=1:size(GMMs.mu,1)
    Pb{1,i}(:,j)= mvnpdf(gesture(:,3:end),GMMs.mu(j,:),GMMs.Sigma(:,:,j));
end

