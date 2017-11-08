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
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Gesture');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Pb');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\States');


%% Emission Probability

% Normalization in order to have the Emission_Prob, as in [6]
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
            Pp = posterior(GMMs{d,1}{f,1},Gesture{d,1}{1,f}(:,4:end));
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

    % 2) Features scaling and mean normalization [5] Grouping [6]
    test=[1,2,4,5];
    
    for c = 1:length(test)
        comp=test(c);

        % 2) Features scaling and mean normalization [5] Grouping [6]
        Total_sequences{d,1} = Scaling_Grouping(Total_sequence{d,1},comp);
        LT=length(Total_sequences{d,1});

        %% Subdivide each user and each trial
        [Data_inusers_rep,LD,K,U,R,T,Gestures]=Subdividing(Total_sequences{d,1});

        %% Beginning Probabilty
        [Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,T,Gestures);

         %% Transition Probabilty
        [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,T,Gestures);

        %% Save

        HMM_tot_generic{d,1}{comp,1}.Beginning_P=Beginning_P;
        HMM_tot_generic{d,1}{comp,1}.Trans_P=Trans_P;
        clear Total_sequences;
    end
end
save test test
save HMM_tot_generic HMM_tot_generic
