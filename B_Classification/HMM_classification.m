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
% 1) Test loading
% [DATA,F,task,user] = Loader_Gestures();
load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\DATA_SETs.mat');
F = 38;
task = 3;

% 2) Loading of the models
% tmp = 'do you want to use a generic model or A specific one? '
% choice = input(tmp);
choice=2;

if choice == 1
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMMs_gestures\HMMs.mat');
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite_generic\HMM_tot_generic.mat');
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite_generic\test.mat');
    
    user = [];
    for d=1:length(HMMs)

        % Initializations
        Accuracy=[];
        Scaling_factor=[];
        Accuracy_mean=[];
        Accuracy_tot=[];
        Accuracy_mean_red=[];
        SUMMARY{d,1}=[];
        user = [user,DATA_SETs{d,1}.test_user];
        DATA = DATA_SETs{d,1}.Test;

        for c = 1:length(test)
            comp=test(c);
            HMM=HMM_tot_generic{d,1}{comp,1};

            if length(HMM)~=0
                for k=1:length(test)
                    SET=test(k);
                    for T=1:task
                        t= ['t',int2str(T)];
                        u = ['u',int2str(length(user))];

                            for i_session=1:numel(fieldnames(DATA.(t)))
                                s = ['s',int2str(i_session)];

                                if F ~=38
                                    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\WLDA');
                                    DATA_reducted.(u).(t).(s) = DATA.(t).(s)(:,3:end-1)*WLDA;
                                    DATA_fin.(u).(t).(s) = [DATA.(t).(s)(:,1:2),DATA_reducted.(u).(t).(s),DATA.(t).(s)(:,end)];
                                else
                                    DATA_fin.(u).(t).(s) = DATA.(t).(s); 
                                end

                                % 3) Features scaling and mean normalization [5] Grouping [6]
                                [New_samples{comp,1}{SET,1}.(u).(t).(s)] = Scaling_Grouping(DATA_fin.(u).(t).(s),SET);
                                LT{comp,1}{SET,1}.(u).(t).(s)=size(New_samples{comp,1}{SET,1}.(u).(t).(s),1);
                                Gesture=size(HMMs{d,1},2);
                                K=3; %num of gaussians in each gesture

                                % 4) CLASSIFICATION
                                Unknown_sequence = New_samples{comp,1}{SET,1}.(u).(t).(s)(:,4:end-1);
                                MultiFactor = HMM.Beginning_P;
                                V=0;
                                Prob=0;

                                for i=1:LT{comp,1}{SET,1}.(u).(t).(s)

                                    Data=Unknown_sequence(i,:);
                                %     [V,Prob,MultiFactor]=Viterbi_Composite(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);
                                %     [V,Prob]=Viterbi_Composite2(HMM,HMMs,Data,i,Gesture,K,V,Prob); % no MultiFactor  
                                %     [V,Prob,MultiFactor]=Viterbi_Composite3(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);  % log % normalizzando Pb 
                                %     [V,Prob,MultiFactor]=Viterbi_Composite4(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob,Sum_tot1); % normalizzando Pb
                                    [V,Prob]=Viterbi_Composite5(HMM,HMMs{d,1},Data,i,Gesture,K,V,Prob); % no MultiFactor  %log

                                    DECODED_SEQUENCE{comp,1}{SET,1}.(u).(t).(s)(i,1)=V;

                                    if Prob == inf || Prob == 0      
                                        break
                                    end

                                end
                                true_T=length(find((DECODED_SEQUENCE{comp,1}{SET,1}.(u).(t).(s)-New_samples{comp,1}{SET,1}.(u).(t).(s)(1:i,end))==0));
                                Accuracy=[Accuracy, true_T/i];
                            end
                    end

                    Accuracy_mean=[ Accuracy_mean; mean(Accuracy)];
                    Accuracy_mean_red=[Accuracy_mean_red;mean(Accuracy(1,find(Accuracy~=0)))];
                    Accuracy_tot=[Accuracy_tot;Accuracy];
                    Accuracy=[];
                    Scaling_factor=[Scaling_factor; SET];
                end
    %                 Accuracy_mean=[ Accuracy_mean; mean(Accuracy)];
    %                 Accuracy_mean_red=[Accuracy_mean_red;mean(Accuracy(1,find(Accuracy~=0)))];
    %                 Accuracy_tot=[Accuracy_tot;Accuracy];
    %                 Accuracy=[];
    %                 Scaling_factor=[Scaling_factor; SET];
            end

                SUMMARY{d,1}=[SUMMARY{d,1}; comp*ones(size(Scaling_factor,1),1),Scaling_factor,Accuracy_tot,Accuracy_mean,Accuracy_mean_red];
                Accuracy=[];
                Accuracy_mean=[];
                Accuracy_tot=[];
                Scaling_factor=[];
                Accuracy_mean_red=[];
        end
    end

    save SUMMARY_generic SUMMARY
else
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMMs_gestures\HMMs.mat');
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite_specific\HMM_tot_specific.mat');
    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\HMM_composite_specific\test.mat');

    user = [];
    for d=1:length(HMMs)

        % Initializations
        Accuracy=[];
        Scaling_factor=[];
        Accuracy_mean=[];
        Accuracy_tot=[];
        Accuracy_mean_red=[];
        user = [user,DATA_SETs{d,1}.test_user];
        DATA = DATA_SETs{d,1}.Test;
        for T=1:task
            SUMMARY{d,1}{T,1}=[];
            for c = 1:length(test)
                comp=test(c);
                    for k=1:length(test)
                        SET=test(k);
                    
                        t= ['t',int2str(T)];
                        u = ['u',int2str(length(user))];
                        HMM=HMM_tot_specific{T,1}{d,1}{comp,1};
                        
                        
                        if length(HMM)~=0

                            for i_session=1:numel(fieldnames(DATA.(t)))
                                s = ['s',int2str(i_session)];

                                if F ~=38
                                    load('C:\Users\Francesco-Greg\Desktop\GMM+HMM\GMM\WLDA');
                                    DATA_reducted.(u).(t).(s) = DATA.(t).(s)(:,3:end-1)*WLDA;
                                    DATA_fin.(u).(t).(s) = [DATA.(t).(s)(:,1:2),DATA_reducted.(u).(t).(s),DATA.(t).(s)(:,end)];
                                else
                                    DATA_fin.(u).(t).(s) = DATA.(t).(s); 
                                end

                                % 3) Features scaling and mean normalization [5] Grouping [6]
                                [New_samples{comp,1}{SET,1}.(u).(t).(s)] = Scaling_Grouping(DATA_fin.(u).(t).(s),SET);
                                LT{comp,1}{SET,1}.(u).(t).(s)=size(New_samples{comp,1}{SET,1}.(u).(t).(s),1);
                                Gesture=size(HMMs{d,1},2);
                                K=3; %num of gaussians in each gesture

                                % 4) CLASSIFICATION
                                Unknown_sequence = New_samples{comp,1}{SET,1}.(u).(t).(s)(:,4:end-1);
                                MultiFactor = HMM.Beginning_P;
                                V=0;
                                Prob=0;

                                for i=1:LT{comp,1}{SET,1}.(u).(t).(s)

                                    Data=Unknown_sequence(i,:);
                                %     [V,Prob,MultiFactor]=Viterbi_Composite(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);
                                %     [V,Prob]=Viterbi_Composite2(HMM,HMMs,Data,i,Gesture,K,V,Prob); % no MultiFactor  
                                %     [V,Prob,MultiFactor]=Viterbi_Composite3(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);  % log % normalizzando Pb 
                                %     [V,Prob,MultiFactor]=Viterbi_Composite4(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob,Sum_tot1); % normalizzando Pb
                                    [V,Prob]=Viterbi_Composite5(HMM,HMMs{d,1},Data,i,Gesture,K,V,Prob); % no MultiFactor  %log

                                    DECODED_SEQUENCE{comp,1}{SET,1}.(u).(t).(s)(i,1)=V;

                                    if Prob == inf || Prob == 0      
                                        break
                                    end

                                end
                                true_T=length(find((DECODED_SEQUENCE{comp,1}{SET,1}.(u).(t).(s)-New_samples{comp,1}{SET,1}.(u).(t).(s)(1:i,end))==0));
                                Accuracy=[Accuracy, true_T/i];
                            end
                        end

                        Accuracy_mean=[ Accuracy_mean; mean(Accuracy)];
                        Accuracy_mean_red=[Accuracy_mean_red;mean(Accuracy(1,find(Accuracy~=0)))];
                        Accuracy_tot=[Accuracy_tot;Accuracy];
                        Accuracy=[];
                        Scaling_factor=[Scaling_factor; SET];
                    end
                    SUMMARY{d,1}{T,1}=[SUMMARY{d,1}{T,1}; comp*ones(size(Scaling_factor,1),1),Scaling_factor,Accuracy_tot,Accuracy_mean,Accuracy_mean_red];
                    Accuracy=[];
                    Accuracy_mean=[];
                    Accuracy_tot=[];
                    Accuracy_mean_red=[];
                    Scaling_factor=[];
    %                 Accuracy_mean=[ Accuracy_mean; mean(Accuracy)];
    %                 Accuracy_mean_red=[Accuracy_mean_red;mean(Accuracy(1,find(Accuracy~=0)))];
    %                 Accuracy_tot=[Accuracy_tot;Accuracy];
    %                 Accuracy=[];
    %                 Scaling_factor=[Scaling_factor; SET];
            end

%                 SUM{d,1}{T,1}=[SUM{d,1}{T,1}; comp*ones(size(Scaling_factor,1),1),Scaling_factor,Accuracy_tot,Accuracy_mean,Accuracy_mean_red];
%                 Accuracy=[];
%                 Accuracy_mean=[];
%                 Accuracy_tot=[];
%                 Scaling_factor=[];
%                 Accuracy_mean_red=[];
        end
    end

    save SUMMARY_specific SUMMARY
end

%% Classification
% % % definition of a random sample
% % Sample_num=randi([0,LT],1,1);
% % Sample_all=New_samples(Sample_num,:);
% % Sample_test=Sample_all(1,3:end-1);
% Unknown_sequence = New_samples(:,3:end-1);
% % Unknown_sequence(358:380,:)=[]; 
% LT=size(Unknown_sequence,1);
% MultiFactor = HMM.Beginning_P;
% V=0;
% Prob=0;
% for i=1:LT
%     Data=Unknown_sequence(i,:);
% %     [V,Prob,MultiFactor]=Viterbi_Composite(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);
% %     [V,Prob]=Viterbi_Composite2(HMM,HMMs,Data,i,Gesture,K,V,Prob); % no MultiFactor  
% %     [V,Prob,MultiFactor]=Viterbi_Composite3(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob);  % log % normalizzando Pb 
% %     [V,Prob,MultiFactor]=Viterbi_Composite4(HMM,HMMs,Data,MultiFactor,i,Gesture,K,V,Prob,Sum_tot1); % normalizzando Pb
%     [V,Prob]=Viterbi_Composite5(HMM,HMMs,Data,i,Gesture,K,V,Prob); % no MultiFactor  %log
%     
%     DECODED_SEQUENCE(i,1)=V;
%      if i== 841
%          l=0;
%      end
% %     figure(1)
% %     plot(P)
%     if Prob == inf || Prob == 0      
%         break
%     end
% 
% end

%% Evaluation
% true_T=length(find((DECODED_SEQUENCE-New_samples(1:i,end))==0));
% Accuracy=true_T/i;

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
