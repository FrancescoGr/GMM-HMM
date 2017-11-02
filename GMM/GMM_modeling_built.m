function [obj]= GMM_modeling_built(Training_set)

%%  Initialization

% 1) Number of Gaussians in the model
K = 3;
obj.AIC = 100;
AIC_old = inf;
d=1;
Iter=1;
Priors_tot=[];

 while ( K == 3 )
     
 % Control variables 
     if K ~= 2
         AIC_old(d,1) = obj.AIC;
         d=d+1;
     end
     check =1;
     CONV=0;

     
     while CONV < 1
        % 1) Computation
        options = statset('Display','final');
        obj = fitgmdist(Training_set,K,'Options',options,'Regularize',0); %Pb
        
        % 2) Definition of a set of solutions:stack "check"; Once, 2 of them are
        % close, one solution is found.
        Covariances_test{1,check} =obj.Sigma;
        Means_test{1,check} = obj.mu;
        Priors_test{1,check} = obj.ComponentProportion;
        [error_Mconv,error_Cconv,error_Pconv]=Err_conv(Means_test,check,Covariances_test,Priors_test);
        
        % check if The best possible solution is "close" enough to the other
        % contained in the stack "check"
        if abs(error_Mconv) < 0.01 && abs(error_Cconv) < 0.01 && abs(error_Pconv) < 0.01
            CONV=1;
        else
            check=check+1;
        end
        
     end
     
     % 3) Stability control: evalutation of the Solution, is it stable
     % enough considering different good solution?
     [OK,Priors_tot]=StabilityControl(Priors_tot,obj.ComponentProportion,Iter);  

    % End of the loop, the stable solution is achived
    if OK 
        K = K+1;
    else
        Iter=Iter+1; %if not, try with a new loop
    end

end