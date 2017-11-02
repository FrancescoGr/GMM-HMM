%% GMMS TRAINING PHASE definition of gaussians that are able to characterize each gesture
% REFERENCES
% [1] Simple Methods for Initializing the EM Algorithm for Gaussian Mixture Models
% [2] The Multivariate Gaussian Distribution (?)
% [3] Bayesian information criterion + Akaike (?)
% [4] OBJECTIVE EVALUATION OF LAPAROSCOPIC SURGICAL SKILLS USING HIDDEN MARKOV 
%     MODELS BASED ON HAPTIC INFORMATION AND TOOL/TISSUE INTERACTIONS
% [5] Objective Skill Evaluation for Laparoscopic Training Based on Motion Analysis
% [6] A Dataset and Benchmarks for Segmentation and Recognition of Gestures in Robotic Surgery

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             preThesis: Francesco Grigoli AA 2016/2017                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear 
close all
clc

% 1) Loading
Data_withGestures = Loader_Gestures();
[r,~]=find(isnan(Data_withGestures));
Data_withGestures(r,:)=[];

Data_withGestures(:,3:end-1)=Data_withGestures(:,3:end-1);
K=max(Data_withGestures(:,end));

% 2) Features scaling and mean normalization [5] Grouping [6]
Data_withGestures = Scaling_Grouping(Data_withGestures);
LT=length(Data_withGestures);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 3) Feature Reduction
% %    - 1st method: PCA
% Training_set = Data_withGestures(:,2:end-1);
% Cross=(1/(size(Training_set,1)))*Training_set'*Training_set;
%  [U,S,V]=svd(Cross);
%  
%  f=1;
%  totS=sum(diag(S));
%  sumS=0;
%  
%  % Total Variance in the remaining Features
%  while (((1-sumS/totS) > 0.1) && (f ~= size(S,2)+1))  
%      sumS=0;
%      f=f+1;     
%      for h=1:f
%          sumS=sumS+S(h,h);
%      end
%  end
%  
%  if f==0
%      fprinf('error');
%  end
%  
% U_red=U(:,1:f);
% 
% % features final
% for i=1:size(Training_set,1)  
%     features_fin(i,:)=U_red'*Training_set(i,:)';   
% end
% 
% % Original reducted Data_set
% Data_withGestures=[Data_withGestures(:,1),features_fin,Data_withGestures(:,end)];
% 
%    - 2st method: LDA as in [6]
% [~, WLDA, ~, ~]=mylda(Data_withGestures(:,3:end-1),Data_withGestures(:,end),9);
% Data_withGestures_reducted = Data_withGestures(:,3:end-1)*WLDA;
% Data_withGestures = [Data_withGestures(:,1:2),Data_withGestures_reducted,Data_withGestures(:,end)];
% save WLDA WLDA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4) Grouping of the data, according to the gestures
[Gesture,States] = Grouping(Data_withGestures,K);

% 5) Definition of States that characterize each Gesture 
%    -> each GMM represent a state in the HMM that characterizes the Gesture
 for i=8 %1: length(Gesture)
     gesture=Gesture{1,i};
     if length(gesture)>200
         
            % built-in function
%               [GMMs] = GMM_modeling_built(gesture(:,3:end));
              [GMMs{i,1}]= GMM_modeling_built_new(gesture(:,3:end));
% 
%               Prior{1,i} = GMMs.ComponentProportion;
%               Covariances{1,i} = GMMs.Sigma;
%               Mean{1,i} = GMMs.mu;
%               [Pb{1,i},Pp{1,i}] = Belonging_Prob( Prior{1,i},Covariances{1,i},Mean{1,i},gesture(:,3:end));
              
              for j=1:size(GMMs{i,1}.mu,1)
                 Pb{1,i}(:,j)= mvnpdf(gesture(:,3:end),GMMs{i,1}.mu(j,:),GMMs{i,1}.Sigma(:,:,j));
              end
              
             % our Function
%          [GMMs(1,i),K(1,i),Prob] = GMM_modeling(gesture(:,2:end));
%          Gesture{1,i}=[Gesture{1,i},Prob.*(GMMs(1,i).Prior)'];
     end
     end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  % CHECK Prior analysis
%  for i=1:length(Prior)
%      if length(Prior{1,i})==0;
%      else
%         P(i,:)=sort(Prior{1,i});
%      end
%  end
%  n=1;
%  for i=1:length(P)
%     if P(i,1)== 0
%         P(n,:)=[]; 
%         n=n+1;
%     end
%  end
%  S=std(P);
%  M=mean(P);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6) Post processing-> over Pb  // NOT WORKING
%    - Standardization
% for f=1:length(Pb)
%     P=Pb{1,f}.*Prior{1,f};
%         for j=1:size(P,2)
%             mu(j)=(1/(size(P,1)))*sum(P(:,j)); % mean of the j-feature
%             standd(j)=std(P(:,j));% std of the j-feature
%         end
% 
%         for i=1:size(P,1)
%             for j=1:size(P,2)
%                 P(i,j)=(P(i,j)-mu(j))/standd(j); % scaling and mean normalization
%             end
%         end
%     Pb_new{1,f}=P;
% end



%% Save 

% 1) For our Function
% save Gesture Gesture
% save GMMs GMMs
% save K K
% save Prob Prob

% 2) For Built-in function
% save GMMs GMMs
% save Pb Pb
% save Gesture Gesture
% save Prior Prior
% save Covariances Covariances
% save Mean Mean 
% save Prob Prob
% save Pb Pb
% save Pp Pp
% save States States
