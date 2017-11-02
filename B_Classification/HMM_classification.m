%% HMM CLASSIFICATION PHASE 
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

%% Initializations
% 1) New data loading
[New_samples,F] = Loader_Gestures();

% 2) Loading of the models
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMMs_gestures\HMMs.mat');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite\HMM.mat');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite\Sum_tot1.mat');

if F ~=38
    load('C:\Users\Francesco-Greg\Documents\MATLAB\HMM\2_JHUISIDATASET\SUPERVISED_CLASS\A_Training\GMM+HMM\GMM\9Features_Builtin\WLDA');
    %    - 2st method: LDA as in [6]
%     [~, WLDA, ~, ~]=mylda(New_samples(:,3:end-1),New_samples(:,end),9);
    New_samples_reducted = New_samples(:,3:end-1)*WLDA;
    New_samples = [New_samples(:,1:2),New_samples_reducted,New_samples(:,end)];
end

% 3) Features scaling and mean normalization [5] Grouping [6]
New_samples = Scaling_Grouping(New_samples);
LT=size(New_samples,1);
Gesture=size(HMMs,2);
K=3; %num of gaussians in each gesture


%% Classification
% % definition of a random sample
% Sample_num=randi([0,LT],1,1);
% Sample_all=New_samples(Sample_num,:);
% Sample_test=Sample_all(1,3:end-1);
Unknown_sequence = New_samples(:,3:end-1);
% Unknown_sequence(358:380,:)=[]; 
LT=size(Unknown_sequence,1);
MultiFactor = HMM.Beginning_P;
V=0;
Prob=0;
for i=1:LT
    Data=Unknown_sequence(i,:);
%     [V,Prob,MultiFactor]=Viterbi_Composite(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);
%     [V,Prob]=Viterbi_Composite2(HMM,HMMs,Data,i,Gesture,K,V,Prob); % no MultiFactor  
%     [V,Prob,MultiFactor]=Viterbi_Composite3(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);  % log % normalizzando Pb 
%     [V,Prob,MultiFactor]=Viterbi_Composite4(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob,Sum_tot1); % normalizzando Pb
    [V,Prob]=Viterbi_Composite5(HMM,HMMs,Data,i,Gesture,K,V,Prob); % no MultiFactor  %log
    
    DECODED_SEQUENCE(i,1)=V;
     if i== 841
         l=0;
     end
%     figure(1)
%     plot(P)
    if Prob == inf || Prob == 0      
        break
    end

end

%% Evaluation
true_T=length(find((DECODED_SEQUENCE-New_samples(1:i,end))==0));
Accuracy=true_T/i;

% %% HMM CLASSIFICATION PHASE 
% % REFERENCES
% % [1] Simple Methods for Initializing the EM Algorithm for Gaussian Mixture Models
% % [2] The Multivariate Gaussian Distribution (?)
% % [3] Bayesian information criterion + Akaike (?)
% % [4] OBJECTIVE EVALUATION OF LAPAROSCOPIC SURGICAL SKILLS USING HIDDEN MARKOV 
% %     MODELS BASED ON HAPTIC INFORMATION AND TOOL/TISSUE INTERACTIONS
% % [5] Objective Skill Evaluation for Laparoscopic Training Based on Motion Analysis
% % [6] A Dataset and Benchmarks for Segmentation and Recognition of Gestures in Robotic Surgery
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %             preThesis: Francesco Grigoli AA 2016/2017                  %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% clear
% close all
% clc
% 
% %% Initializations
% % 1) New data loading
% [New_samples,F,sessions] = Loader_Gestures();
% 
% % 2) Loading of the models
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMMs_gestures\HMMs.mat');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite\HMM.mat');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite\Sum_tot1.mat');
% 
% if F ~=38
%     load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\9Features_Builtin\WLDA');
%     %    - 2st method: LDA as in [6]
% %     [~, WLDA, ~, ~]=mylda(New_samples(:,3:end-1),New_samples(:,end),9);
%     New_samples_reducted = New_samples(:,3:end-1)*WLDA;
%     New_samples = [New_samples(:,1:2),New_samples_reducted,New_samples(:,end)];
% end
% 
% % 3) Features scaling and mean normalization [5] Grouping [6]
% for i_session = 1:length(sessions)
%     s = ['s',int2str(i_session)];
%     New_samples.(s) = Scaling_Grouping(New_samples.(s));
%     LT.(s)=size(New_samples.(s),1);
% end
% Gesture=size(HMMs,2);
% K=3; %num of gaussians in each gesture
% 
% 
% 
% %% Classification
% % % definition of a random sample
% % Sample_num=randi([0,LT],1,1);
% % Sample_all=New_samples(Sample_num,:);
% % Sample_test=Sample_all(1,3:end-1);
% % Unknown_sequence = New_samples(:,3:end-1);
% % MultiFactor = HMM.Beginning_P;
% % V=0;
% % Prob=0;
% for i_session = 1:length(sessions)
%     s = ['s',int2str(i_session)];
%     Unknown_sequence = New_samples.(s)(:,3:end-1);
%     MultiFactor = HMM.Beginning_P;
%     V=0;
%     Prob=0;
%     for i=1:LT.(s)
%         Data=Unknown_sequence(i,:);
%         [V,Prob,MultiFactor]=Viterbi_Composite(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);
% %         [V,Prob,MultiFactor]=Viterbi_Composite2(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob); % no MultiFactor  
%     %     [V,Prob,MultiFactor]=Viterbi_Composite3(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);  % log
%     %     [V,Prob,MultiFactor]=Viterbi_Composite4(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob,Sum_tot1); % normalizzando Pb
%         DECODED_SEQUENCE.(s)(i,1)=V;
%     end
% end
% 
% 
% %% Evaluation
% Sum=0;
% true_T=0;
% 
% for i_session = 1:length(sessions)
%     s = ['s',int2str(i_session)];
%     good.(s)=length(find((DECODED_SEQUENCE.(s)-New_samples.(s)(:,end))==0));
%     true_T=true_T+length(find((DECODED_SEQUENCE.(s)-New_samples.(s)(:,end))==0));
%     Sum=Sum+LT.(s);
%     true_S.(s)=good.(s)/LT.(s);
% end
% per=true_T/Sum
