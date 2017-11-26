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
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\choice');

num_tasks = 3;

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

if choice == 1
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

        % Data_withGestures(ind,:)=[];
        % Merging the tola sequence with Data_withGestures

        % 2) Features scaling and mean normalization [5] Grouping [6]
        Total_sequence{d,1} = Scaling_Grouping(Total_sequence{d,1});
        LT=length(Total_sequence{d,1});

        %% Subdivide each user and each trial
        [Data_inusers_rep,LD,K,U,R,T,Gestures]=Subdividing(Total_sequence{d,1});

        %% Beginning Probabilty

        [Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,T,Gestures);

         %% Transition Probabilty

        [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,T,Gestures);

        %% Save

        for f=1:length(Gesture{d,1})
            if f ~=7
                HMMs{d,1}(1,f).Beginning_P=Beginning_P{1,f};
                HMMs{d,1}(1,f).Trans_P=Trans_P{1,f};
                HMMs{d,1}(1,f).Emission_P=Emission_P_sequence{d,1}{1,f};
                HMMs{d,1}(1,f).Prior=GMMs{d,1}{f,1}.ComponentProportion;
                HMMs{d,1}(1,f).Covariances=GMMs{d,1}{f,1}.Sigma;
                HMMs{d,1}(1,f).Mean=GMMs{d,1}{f,1}.mu;
            end
        end
    end 
else
%% Emission Probability
    % Normalization in order to have the Emission_Prob, as in [6]
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Gest_task');
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\Gest_task_ind');
    for d=1:length(GMMs)
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
            
            Total_sequence{d,1}{T,1} = [];
            
            for lab=1:length(G)
                f=G(1,lab);
                P=Pb{d,1}{T,1}{1,f}.*GMMs{d,1}{T,1}{f,1}.ComponentProportion;
                Emission_P_sequence{d,1}{T,1}{1,f}=P;
                Pp = posterior(GMMs{d,1}{T,1}{f,1},Gest_task{d,1}{T,1}{1,f}(:,4:end));
                [~,ind]=max(Pp');
                State{d,1}{T,1}{1,f}=ind';
                            
                % total Sequence of the states and of the gestures
                
                lines = States{d,1}{1,f}(Gest_task_ind{d,1}{T,1}{1,f},1);
                Points =[Gest_task{d,1}{T,1}{1,f},f*ones(size(Gest_task{d,1}{T,1}{1,f},1),1),Emission_P_sequence{d,1}{T,1}{1,f},State{d,1}{T,1}{1,f}];
                if size(Points,1)> 200  % elimina gesture 7,10 pochi samples
                    Total_sequence{d,1}{T,1}(lines,:)=Points;   
                end
            end

            % Deleting of the empty samples (gesture 7,10)
            [ind]=find(Total_sequence{d,1}{T,1}(:,1)==0);
            Total_sequence{d,1}{T,1}(ind,:)=[];

            % Data_withGestures(ind,:)=[];
            % Merging the tola sequence with Data_withGestures

            % 2) Features scaling and mean normalization [5] Grouping [6]
            Total_sequence{d,1}{T,1} = Scaling_Grouping(Total_sequence{d,1}{T,1});
            LT=length(Total_sequence{d,1}{T,1});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% Subdivide each user and each trial
            [Data_inusers_rep,LD,K,U,R,T,Gestures]=Subdividing(Total_sequence{d,1}{T,1});

            %% Beginning Probabilty

            [Beginning_P] = Startprob(Data_inusers_rep,LD,K,U,R,T,Gestures);

             %% Transition Probabilty

            [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,T,Gestures);

            %% Save
            for count=1:Gestures %length(G)
                if length(find(G==count))~=0
                    
                    HMMs{d,1}{T,1}(1,count).Beginning_P =Beginning_P{1,count};
                    HMMs{d,1}{T,1}(1,count).Trans_P=Trans_P{1,count};
                    HMMs{d,1}{T,1}(1,count).Emission_P=Emission_P_sequence{d,1}{T,1}{1,count};                
                    HMMs{d,1}{T,1}(1,count).Prior = GMMs{d,1}{T,1}{count,1}.ComponentProportion;
                    HMMs{d,1}{T,1}(1,count).Covariances=GMMs{d,1}{T,1}{count,1}.Sigma;
                    HMMs{d,1}{T,1}(1,count).Mean=GMMs{d,1}{T,1}{count,1}.mu;

                else
                    HMMs{d,1}{T,1}(1,count).Beginning_P =[];
                    HMMs{d,1}{T,1}(1,count).Trans_P=[];
                    HMMs{d,1}{T,1}(1,count).Emission_P=[];
                    HMMs{d,1}{T,1}(1,count).Prior=[];
                    HMMs{d,1}{T,1}(1,count).Covariances=[];
                    HMMs{d,1}{T,1}(1,count).Mean=[];

                end
                    
            end

        end 
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save HMMs HMMs