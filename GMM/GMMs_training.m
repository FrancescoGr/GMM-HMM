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

% % 1) Loading
[Data_inusersANDtasks,users] = Loader_Gestures();

% 2) Test set and Training Set definition
User_out = [1,2,3]; %['D','C','B']; 
DATA_SETs = Dataset(User_out,Data_inusersANDtasks,users);
save DATA_SETs DATA_SETs
% load('DATA_SETs.mat')


% 3) number of features 
dim = 38;

for d=1:length(DATA_SETs)
    Data_withGestures = DATA_SETs{d,1}.Training;

    [r,~]=find(isnan(Data_withGestures));
    Data_withGestures(r,:)=[];

    K=max(Data_withGestures(:,end));

    % 4) Features scaling and mean normalization [5] Grouping [6]
    Data_withGestures = Scaling_Grouping(Data_withGestures);
    LT=length(Data_withGestures);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % 5) Feature Reduction
    % %   - 5.1st method: PCA
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %   - 5.2st method: LDA as in [6]
    if dim ~= 38
        [~, WLDA{d,1}, ~, ~]=mylda(Data_withGestures(:,4:end-1),Data_withGestures(:,end),dim);
        Data_withGestures_reducted = Data_withGestures(:,4:end-1)*WLDA{d,1};
        Data_withGestures = [Data_withGestures(:,1:3),Data_withGestures_reducted,Data_withGestures(:,end)];
    end

    % 6) Grouping of the data, according to the gestures
    [Gesture{d,1},States{d,1}] = Grouping(Data_withGestures,K);

    % 7) Definition of States that characterize each Gesture 
    %    -> each GMM represent a state in the HMM that characterizes the Gesture
     for i=1: length(Gesture{d,1})
         gesture=Gesture{d,1}{1,i};
         if length(gesture)>200

                % built-in function
    %               [GMMs] = GMM_modeling_built(gesture(:,3:end));
    %               Prior{1,i} = GMMs.ComponentProportion;
    %               Covariances{1,i} = GMMs.Sigma;
    %               Mean{1,i} = GMMs.mu;
    %               [Pb{1,i},Pp{1,i}] = Belonging_Prob( Prior{1,i},Covariances{1,i},Mean{1,i},gesture(:,3:end));
                  
                  [GMMs{d,1}{i,1}]= GMM_modeling_built_new(gesture(:,4:end));
                  for j=1:size(GMMs{d,1}{i,1}.mu,1)
                     Pb{d,1}{1,i}(:,j)= mvnpdf(gesture(:,4:end),GMMs{d,1}{i,1}.mu(j,:),GMMs{d,1}{i,1}.Sigma(:,:,j));
                  end

                 % our Function
    %          [GMMs(1,i),K(1,i),Prob] = GMM_modeling(gesture(:,2:end));
    %          Gesture{1,i}=[Gesture{1,i},Prob.*(GMMs(1,i).Prior)'];
         end
     end
end

 save Pb Pb
 save GMMs GMMs
 save Gesture Gesture
 save States States
 if dim ~= 38
    save WLDA WLDA
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
