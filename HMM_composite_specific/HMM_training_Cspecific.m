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
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Gest_task');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Pb');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\States');
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\DATA_SETs');

%% Initializations
num_tasks = 3;


%% Emission Probability

% Normalization in order to have the Emission_Prob, as in [6]
for T= 1:num_tasks
   % definition of the gestures in each task
    if T==1
       G=[1,11:15];
    end
    if T==2
       G=[1:6,8,11];
    end
    if T==3
       G=[1:6,8:11];
    end
    
    for d=1:length(GMMs)
        
        Total_sequence{d,1}=[];
        for lab=1:length(G)
            f=G(1,lab);
            if f ~=7
                [Gest_task_ind] = find(Gesture{d,1}{1,f}(:,2)== T);
                Gest_task{T,1}{d,1}{1,f} = Gesture{d,1}{1,f}(Gest_task_ind,:);
                 
%                 P=Pb{d,1}{1,f}.*GMMs{d,1}{f,1}.ComponentProportion;
%                 Emission_P_sequence{d,1}{1,f}=P;
% 
%                 Pp = posterior(GMMs{d,1}{f,1},Gesture{d,1}{1,f}(:,4:end));
%                 [~,ind]=max(Pp');
%                 State{d,1}{1,f}=ind';

                lines = States{d,1}{1,f}(Gest_task_ind,:);
                Points =[Gest_task{T,1}{d,1}{1,f},f*ones(size(Gest_task{T,1}{d,1}{1,f},1),1)]; %,1),Emission_P_sequence{d,1}{1,f},State{d,1}{1,f}];
                if size(Points,1)> 200  % elimina gesture 10, pochi samples
                    Total_sequence{T,1}{d,1}(lines,:)=Points;   
                end
            end
        end

    % Deleting of the empty samples (gesture 10)
        [ind]=find(Total_sequence{T,1}{d,1}(:,1)==0);
        Total_sequence{T,1}{d,1}(ind,:)=[];

        % 2) Features scaling and mean normalization [5] Grouping [6]
        test=[1,2,4,5];
        
        for c = 1:length(test)
            comp=test(c);

            % 2) Features scaling and mean normalization [5] Grouping [6]
            Total_sequences{d,1} = Scaling_Grouping(Total_sequence{T,1}{d,1},comp);
            LT=length(Total_sequences{d,1});

            %% Subdivide each user and each trial
            [Data_inusers_rep,LD,K,U,R,Y,Gestures]=Subdividing(Total_sequences{d,1});

            %% Beginning Probabilty
            [Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,Y,Gestures);

             %% Transition Probabilty
            [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,Y,Gestures);

            %% Save

            HMM_tot_specific{T,1}{d,1}{comp,1}.Beginning_P=Beginning_P;
            HMM_tot_specific{T,1}{d,1}{comp,1}.Trans_P=Trans_P;
            clear Total_sequences;
        end
        
    end
end
save test test
save HMM_tot_specific HMM_tot_specific
