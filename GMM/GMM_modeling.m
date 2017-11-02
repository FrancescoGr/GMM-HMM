function [GMMs,K,Prob]= GMM_modeling(Training_set)

LT=length(Training_set);

if LT ~= 0
%%  Initialization

% 1) Number of Gaussians in the model
K = 3;
AIC_new = 0;
AIC = inf;

% while ( K == 2 ||  AIC_new < AIC)
while ( K == 3) % Number of gaussian in each HMM [6]
     
 % Control variables 
     if K ~= 2
         AIC = AIC_new;
     end
     check =1;
     iteration =1;
     CONV=0;
     tic;

     
     while CONV < 1 

    % 2) E-M Initialization [1]

        % Partitions definition
        [Clusters] = Partitions(Training_set,K);

        % Gaussians Parameter
        [Covariances,Means,Priors] = Gaussians(K,Training_set,Clusters);

        % As I though
        n = size(Means,2);

    %     % according to the paper [4], it makes no sense
    %     n = size(Training_set,1);
    
             if length(find(isnan(Covariances)))~=0 || length(find((Covariances)==inf)) ~=0 % NaN nella covaricanza
                 CONV=10 ;
                 control = 0;
             end


    % 3) E-M Process, Model definition 

        iter = 1;
        L = 0;
        Confidence = 1000;

    % 4) Iteretions until the model converges

        while (Confidence(end) > 0.1)% && (iter < 500)
    %         time1=tic;

            % Belonging probability GMM [2]
                for j = 1:length(Priors)
                    for i = 1:LT
                        Pb(i,j) = (1/((sqrt((2*pi)^(n)*det(Covariances(:,:,j))))))*exp(-0.5*(Training_set(i,:)-Means(j,:))*pinv(Covariances(:,:,j))*(Training_set(i,:)-Means(j,:))');
                    end

                end

            % Posterior Probability (Baysian)
                for i = 1:LT

                    Sum = 0;
                    for k = 1:K
                        Sum = Sum+(Pb(i,k)*Priors(k,1));
                    end

                    for j = 1:length(Priors)
                        if Sum == 0
                            l=0;
                        end
                        Pp(i,j) = (Pb(i,j)*Priors(j,1))/Sum;
                    end     
                end

            % Weights definition (necessary to update the parameters)
              for  j = 1:length(Priors)
                  for i = 1:LT
                      Weights(i,j) = Pp(i,j)/sum(Pp(:,j));
                  end   
              end

            % Gaussians Parameters updating
            Covariances_new = zeros(size(Training_set,2),size(Training_set,2),K);

             for  j = 1:length(Priors)

                 % Priors
                 Priors_new(j,1) = (1/length(Training_set))*sum(Pp(:,j));   

                 % Means 
                 for i = 1:size(Means,2)
                     Means_new(i,j) = Weights(:,j)'*Training_set(:,i) ;
                 end

                 % Gaussians
                Covariances_new = Covariance(Training_set,Means_new,Covariances_new,j,Weights);

             end

             % Convergence evaluation
             L_new = Evaluation(LT,Pb,Priors);
             Confidence(iter,1) = abs (L_new-L);
             
%              if iter>2   %local minimum
%                 if  Confidence(end,1)== Confidence(end-1,1)
% 
%                     %  E-M Initialization
%                     % Partitions definition
%                     [Clusters] = Partitions(Training_set,K);
% 
%                     % Gaussians Parameter
%                     [Covariances,Means,Priors] = Gaussians(K,Training_set,Clusters);
% 
%                     % E-M Process, Model definition 
%                     iter = 1;
%                     L = 0;
%                     Confidence = 1000;
%                     fprintf('restart 1\n')
%                   end
%               end


             if length(find(isnan(Covariances_new)))~=0 % NaN nella covaricanza
                 iter=401 ;
                 Confidence(end,1)= 121;
             end

    %          time1=toc;
    %          fprintf('iter= %d time=%f Confidence=%f \n',iter,time1,Confidence(end,1))
             fprintf('iter= %d Confidence=%f \n',iter,Confidence(end,1))


             Priors = Priors_new;
             Means = Means_new';
             Covariances = Covariances_new;
             L = L_new;

             iter=iter+1;
             iteration =iteration+1;
             
              if size(Confidence,1)>2   %local minimum
                  if  Confidence(end,1) == Confidence(end-1,1)

                    %  E-M Initialization
                    % Partitions definition
                    [Clusters] = Partitions(Training_set,K);

                    % Gaussians Parameter
                    [Covariances,Means,Priors] = Gaussians(K,Training_set,Clusters);

                    % E-M Process, Model definition 
                    iter = 1;
                    L = 0;
                    Confidence = 1000;
                    fprintf('restart 1\n')
                  end
              end
                  
             if iter>100 && Confidence(end,1)> 60 
                %  E-M Initialization
                % Partitions definition
                [Clusters] = Partitions(Training_set,K);

                % Gaussians Parameter
                [Covariances,Means,Priors] = Gaussians(K,Training_set,Clusters);

                % E-M Process, Model definition 
                iter = 1;
                L = 0;
                Confidence = 1000;
                fprintf('restart 2\n')
             end
             control = 1;
        end
        
        
% Check the Convergence
        if control == 1
             Covariances_test{1,check} =Covariances;
             Means_test{1,check} = Means;
             Priors_test{1,check} = Priors;
             [error_Mconv,error_Cconv,error_Pconv]=Err_conv(Means_test,check,Covariances_test,Priors_test);

        if abs(error_Mconv) < 0.01 && abs(error_Cconv) < 0.1 && error_Pconv < 0.01
                CONV=1;
             else
                check=check+1;
             end
        end
     end
     
    clear   Covariances_test ;
    clear   Means_test;
% 5) Number of Gaussians evaluation [3]
    
    % Number of parameters in the model
    DIM = length(Covariances);
    p = K*((DIM*DIM-DIM)/2+2*DIM+1)-1;              % p OK

    % Akaike Criterion
    [AIC_new,BIC_new] = Akaike(L,p,LT);  % BIC_new OK
    time=toc;
%     
% % 6) Results
%      [Accuracy_def{1,K},index_class1] = test(Training_set,LT,Covariances,Means,Priors,Pp,K);
%      fprintf('Classification results: EM_ALGORITHM\n');
%      fprintf(' Accuracy = %f error_conv=%f \n',Accuracy_def{1,K},error_Mconv);
%      fprintf('%f iterations ,log-likelihood = %f \n',iteration,L);
%      fprintf('BIC = %f ,AIC = %f,time=%f \n\n',BIC_new,AIC_new,time);
%      
%      x1 = -6:1:6; x2= -6:1:6;
%      [X1,X2] = meshgrid(x1,x2);
%  
%      hold on        
%      for i=1:K
%         F = mvnpdf([X1(:) X2(:)],Means(i,:),Covariances(:,:,i));
%         F = reshape(F,length(x2),length(x1));
%         c=contour(x1,x2,F);
%         hold on
%      end

          %% CLASSIFICATION WITH THE BUILT-IN FUNCTION

% % 1) Computation
% tic;
% options = statset('Display','final');
% fprintf('Classification results: built-in function \n');
% obj = fitgmdist(Training_set,K,'Options',options);
% time =toc;
% 
% % 2) Results
%  [Accuracy,index_class2]= test(Training_set,LT,obj.Sigma,obj.mu,obj.ComponentProportion',Pp,K);
%  [Classification_error,non_errors,n1,n2]=ClassError(index_class1,index_class2);
% %  fprintf(' Accuracy = %f \n',Accuracy); % not usable, computed over the  prediction made by EM
%  fprintf('BIC = %f ,AIC = %f,time=%f \n\n',obj.BIC,obj.AIC,time);
%  fprintf('TOTAL ERROR (emVSbuilt-in)=%f, non_errors=%d over %d trans_EM and %d trans_test\n\n',Classification_error,non_errors,n1,n2);

 % figure();
%  scatter(Training_set(:,1),Training_set(:,2),10,'.');
%  hold on
%  h = ezcontour(@(x,y)pdf(obj,[x y]),[-6 6],[-6 6]);

     
    
         %% Saving the Variables    

Priors_def{1,K}=Priors;
Means_def{1,K}=Means;
Covariances_def {1,K}= Covariances;
% L_def{1,K}= L;
Pb_def{1,K}=Pb ;
% Pp_def{1,K}=Pp ;
% AIC_def{1,K} = AIC_new;
% BIC_def{1,K} = BIC_new;
% Confidence_def{1,K}=Confidence;
% 
% save Priors_def Priors_def;
% save Means_def Means_def;
% save Covariances_def  Covariances_def ;
% save L_def L_def;
% save Pb_def Pb_def;
% save Pp_def Pp_def;
% save AIC_def AIC_def;
% save BIC_def BIC_def;
% save Confidence_def Confidence_def;
% save Accuracy_def Accuracy_def;



%     new loop if the model it has been improved
%     if ( BIC_new > BIC || AIC_new < AIC)
    if ( AIC_new < AIC || K==2) 
        K = K+1
    end
    
end
end

K=K-1;

GMMs.Prior = Priors_def{1,K};
GMMs.Means=Means_def{1,K};
GMMs.Covariances=Covariances_def{1,K};
Prob = Pb_def{1,K};

end