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

% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Mean');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Pp');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prior');
% load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\AllFeatures_Builtin\Prob');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Pb');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\GMMs');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Gesture');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\States');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\DATA_SETs');


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
for d=1:length(GMMs)

    for f=1:length(Gesture{d,1})
        if f ~=7
            P=Pb{d,1}{1,f}.*GMMs{d,1}{f,1}.ComponentProportion;
        %     tot=sum(P);
        %     P=P./tot;
            Emission_P_sequence{d,1}{1,f}=P;
        end
    %     Sum_tot{1,f}(:,:)=tot;
    end

%% Sequence reconstruction

    for f=1:length(Gesture{d,1})
        if f ~=7
            Pp = posterior(GMMs{d,1}{f,1},Gesture{d,1}{1,f}(:,3:end));
            [~,ind]=max(Pp');
            State{d,1}{1,f}=ind';
        end
    end

    % total Sequence of the states and of the gestures
    Total_sequence{d,1}=[];
    for f=1:length(Gesture{d,1})
        if f ~=7
            lines = States{d,1}{1,f};
            Points =[Gesture{d,1}{1,f},f*ones(size(Gesture{d,1}{1,f},1),1),Emission_P_sequence{d,1}{1,f},State{d,1}{1,f}];
            if size(Points,1)> 200  % elimina gesture 7,10 pochi samples
                Total_sequence{d,1}(lines,:)=Points;   
            end
        end
    end

    % Deleting of the empty samples (gesture 7,10)
    [ind]=find(Total_sequence{d,1}(:,1)==0);
    Total_sequence{d,1}(ind,:)=[];

    % Data_withGestures(ind,:)=[];
    % Merging the tola sequence with Data_withGestures

    % 2) Features scaling and mean normalization [5] Grouping [6]
    Total_sequence{d,1} = Scaling_Grouping(Total_sequence{d,1});
    LT=length(Total_sequence{d,1});

    %% Subdivide each user and each trial
    [Data_inusers_rep,LD,K,U,R,Gestures]=Subdividing(Total_sequence{d,1});

    %% Beginning Probabilty

    [Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,Gestures);

     %% Transition Probabilty

    [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,Gestures);

    %% Save

    for f=1:length(Gesture)
        if f ~=7
            HMMs{d,1}(1,f).Beginning_P=Beginning_P{1,f};
            HMMs{d,1}(1,f).Trans_P=Trans_P{1,f};
            HMMs{d,1}(1,f).Emission_P=Emission_P_sequence{1,f};
            HMMs{d,1}(1,f).Prior=GMMs{f,1}.ComponentProportion;
            HMMs{d,1}(1,f).Covariances=GMMs{f,1}.Sigma;
            HMMs{d,1}(1,f).Mean=GMMs{f,1}.mu;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save HMMs HMMs