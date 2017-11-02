%% HMMs TRAINING PHASE 
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

load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\GMMs');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Gesture');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Mean');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Pb');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Pp');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prior');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prob');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\States');

% % 1) Loading
% Data_withGestures = Loader_Gestures();
% [r,~]=find(isnan(Data_withGestures));
% Data_withGestures(r,:)=[];
% 
% Data_withGestures(:,3:end-1)=Data_withGestures(:,3:end-1);
% K=max(Data_withGestures(:,end));
% 
% % 2) Features scaling and mean normalization [5] Grouping [6]
% Data_withGestures = Scaling_Grouping(Data_withGestures);
% LT=length(Data_withGestures);

%% Emission Probability

% Normalization in order to have the Emission_Prob, as in [6]
for f=1:length(Gesture)
    P=Pb{1,f}.*GMMs{f,1}.ComponentProportion;
%     tot=sum(P);
%     P=P./tot;
    Emission_P_sequence{1,f}=P;
%     Sum_tot{1,f}(:,:)=tot;
end

%% Sequence reconstruction

for f=1:length(Gesture)
    Pp = posterior(GMMs{f,1},Gesture{1,f});
    [~,ind]=max(Pp');
    State{1,f}=ind';
end

% total Sequence of the states and of the gestures
Total_sequence=[];
for f=1:length(Gesture)
    lines = States{1,f};
    Points =[Gesture{1,f},f*ones(size(Gesture{1,f},1),1),Emission_P_sequence{1,f},State{1,f}];
    if size(Points,1)> 200  % elimina gesture 10, pochi samples
        Total_sequence(lines,:)=Points;   
    end
end

% Deleting of the empty samples (gesture 10)
[ind]=find(Total_sequence(:,1)==0);
Total_sequence(ind,:)=[];

% Data_withGestures(ind,:)=[];
% Merging the tola sequence with Data_withGestures

% 2) Features scaling and mean normalization [5] Grouping [6]
Total_sequence = Scaling_Grouping(Total_sequence);
LT=length(Total_sequence);

%% Subdivide each user and each trial
[Data_inusers_rep,LD,K,U,R,Gestures]=Subdividing(Total_sequence);

%% Beginning Probabilty

[Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,Gestures);

 %% Transition Probabilty

[Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,Gestures);

%% Save

for f=1:length(Gesture)
    HMMs(1,f).Beginning_P=Beginning_P{1,f};
    HMMs(1,f).Trans_P=Trans_P{1,f};
    HMMs(1,f).Emission_P=Emission_P_sequence{1,f};
    HMMs(1,f).Prior=GMMs{f,1}.ComponentProportion;
    HMMs(1,f).Covariances=GMMs{f,1}.Sigma;
    HMMs(1,f).Mean=GMMs{f,1}.mu;
end


save HMMs HMMs