%% COMPOSITE HMM, training phase
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

% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Covariances');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Gesture');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Mean');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Pb');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Pp');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prior');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prob');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\States');

load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\GMMs');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Gesture');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Mean');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Pb');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Pp');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prior');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prob');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\States');


%% Emission Probability

% Normalization in order to have the Emission_Prob, as in [6]
Sum_tot1=zeros(1,3);
for f=1:length(Gesture)
    if f ~=7
        P=Pb{1,f}.*GMMs{f,1}.ComponentProportion;
        tot=sum(P);
        P=P./tot;
        Emission_P_sequence{1,f}=P;
        Sum_tot{1,f}(:,:)=tot;
        Sum_tot1=Sum_tot1+tot;
    end
end
Sum_tot1=Sum_tot1./length(Gesture);

save Sum_tot1 Sum_tot1

%% Sequence reconstruction

% for f=1:length(Gesture)
%     [~,ind]=max(Pp{1,f}');
%     State{1,f}=ind';
% end


for f=1:length(Gesture)
    if f ~=7
        Pp = posterior(GMMs{f,1},Gesture{1,f}(:,3:end));
        [~,ind]=max(Pp');
        State{1,f}=ind';
    end
end


% total Sequence
Total_sequence=[];
for f=1:length(Gesture)
    if f ~=7
        lines = States{1,f};
        Points =[Gesture{1,f},f*ones(size(Gesture{1,f},1),1),Emission_P_sequence{1,f},State{1,f}];
        if size(Points,1)> 200  % elimina gesture 10, pochi samples
            Total_sequences(lines,:)=Points;   
        end
    end
end

% Deleting of the empty samples (gesture 10)
[ind]=find(Total_sequences(:,1)==0);
Total_sequences(ind,:)=[];

% 2) Features scaling and mean normalization [5] Grouping [6]
test=[1,2,4,5];
for c = 1:length(test)
    comp=test(c);

    Total_sequence = Scaling_Grouping(Total_sequences,comp);
    LT=length(Total_sequence);

    %% Subdivide each user and each trial
    [Data_inusers_rep,LD,K,U,R,Gestures]=Subdividing(Total_sequence);

    %% Beginning Probabilty

    [Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,Gestures);

     %% Transition Probabilty
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,Gestures);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Save

    HMM_tot{comp,1}.Beginning_P=Beginning_P;
    HMM_tot{comp,1}.Trans_P=Trans_P;
    clear Total_sequence;
end

save test test
save HMM_tot HMM_tot
