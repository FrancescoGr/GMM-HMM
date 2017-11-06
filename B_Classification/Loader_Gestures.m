function [Features,F,task,user] = Loader_Gestures()


task=3;
user='F';

% tmp = 'feataures under consideration = ';
F =38; %input(tmp);

%% Gestures
%      sessions= [2];
%      user = ['E']; 
for T=1:task
    if T==1
        dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set/Knot_Tying_';
        dir2='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set\transcriptions/Knot_Tying_';
    end
    if T==2
        dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set/Needle_Passing_';
        dir2='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set\transcriptions/Needle_Passing_';
    end
    if T==3     
        dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set/Suturing_';
        dir2='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set\transcriptions/Suturing_';
    end

    for j=1:length(user)
        u = ['u',int2str(j)];
        if T==1 || T==3
            sessions=[1:5];
        else
            sessions=[1,3,4]; 
        end
        for i=1:length(sessions)
            i_session= sessions(i);
            s = ['s',int2str(i_session)];
            filename = [dir2,user,num2str(0),num2str(0),num2str(i_session),'.txt'];
            raw_data.(u).(s) = ExtractData_gestures(filename);
        end
    end


    %% Pairing Gestures to real data

    for j=1:length(user)
        u = ['u',int2str(j)];
        if T==1 || T==3
            sessions=[1:5];
        else
            sessions=[1,3,4]; 
        end
        for i=1:length(sessions)
            i_session= sessions(i);
            s = ['s',int2str(i_session)];
            filename = [dir1,user,num2str(0),num2str(0),num2str(i_session),'.txt'];
            data.(u).(s) = ExtractData(filename);
        end
    end
        


    % Features extraction (useful attributes)

    for j=1:length(user)
         u = ['u',int2str(j)];
         if T==1 || T==3
             sessions=[1:5];
         else
             sessions=[1,3,4]; 
         end
         t= ['t',int2str(T)];
         for i=1:length(sessions)
             i_session= sessions(i);
             s = ['s',int2str(i_session)];
             Features.(u).(t).(s)=[];
             Session=[];
             Session=[Session,j*(ones(length(data.(u).(s).PSM1.pos),1))];
             Session=[Session,i_session*(ones(length(data.(u).(s).PSM1.pos),1))];     
             Session=[Session,data.(u).(s).PSM1.pos];
             Session=[Session,data.(u).(s).PSM1.rot]; 
             Session=[Session,data.(u).(s).PSM1.vel]; 
             Session=[Session,data.(u).(s).PSM1.rot_vel]; 
             Session=[Session,data.(u).(s).PSM1.grip_vel]; 
             Session=[Session,data.(u).(s).PSM2.pos];
             Session=[Session,data.(u).(s).PSM2.rot]; 
             Session=[Session,data.(u).(s).PSM2.vel]; 
             Session=[Session,data.(u).(s).PSM2.rot_vel]; 
             Session=[Session,data.(u).(s).PSM2.grip_vel]; 

             k=size(Session,2);
        % Adding the gesture
             for m=1:length(raw_data.(u).(s).start)
                 b=raw_data.(u).(s).start(m,1);
                 e=raw_data.(u).(s).stop(m,1);
                 Session(b:e,k+1)=raw_data.(u).(s).gestures(m,1);
             end
        % Cut of the session  according to its gesture
             Session=Session(raw_data.(u).(s).start(1):raw_data.(u).(s).stop(end),:);            
             if length(find(Session(:,end)==0))~=0
                 ind=find(Session(:,end)==0);
                 Session(ind,:)=[];
             end
                
             Features.(u).(t).(s) = [Session];

         end 
    
    end
end
end
% function [Features,F] = Loader_Gestures()


% tmp = 'Kind of tasks = ';
% T = input(tmp);
% 
% tmp = 'User = ';
% user = input(tmp,'s');
% 
% tmp = 'Session = ';
% sessions = input(tmp);
% T=[1,2,3]
% user='F';
% sessions =[1:5];
% 
% % tmp = 'feataures under consideration = ';
% F =38; %input(tmp);

%% Gestures
%      sessions= [2];
%      user = ['E']; 

% if T==1
%     dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set/Knot_Tying_';
%     dir2='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set\transcriptions/Knot_Tying_';
% end
% if T==2
%     dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set/Needle_Passing_';
%     dir2='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set\transcriptions/Needle_Passing_';
% end
% if T==3     
%     dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set/Suturing_';
%     dir2='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Test_Set\transcriptions/Suturing_';
% end
% 
%     for j=1:length(user)
%         u = ['u',int2str(j)];
%         for i_session = sessions
%             s = ['s',int2str(i_session)];
%             filename = [dir2,user,num2str(0),num2str(0),num2str(i_session),'.txt'];
%             raw_data.(u).(s) = ExtractData_gestures(filename);
%         end
%     end
% 
% 
%     %% Pairing Gestures to real data
%         for j=1:length(user)
%         u = ['u',int2str(j)];
%             for i_session = sessions
%                 s = ['s',int2str(i_session)];
%                 filename = [dir1,user,num2str(0),num2str(0),num2str(i_session),'.txt'];
%                 data.(u).(s) = ExtractData(filename);
%             end
%         end
%         
% 
% 
%     % Features extraction (useful attributes)
%     Features=[];
% 
%         for j=1:length(user)
%             u = ['u',int2str(j)];
%             for i_session = sessions
% 
%                 s = ['s',int2str(i_session)];
%                 Session=[];
%                 Session=[Session,j*(ones(length(data.(u).(s).PSM1.pos),1))];
%                 Session=[Session,i_session*(ones(length(data.(u).(s).PSM1.pos),1))];     
%                 Session=[Session,data.(u).(s).PSM1.pos];
%                 Session=[Session,data.(u).(s).PSM1.rot]; 
%                 Session=[Session,data.(u).(s).PSM1.vel]; 
%                 Session=[Session,data.(u).(s).PSM1.rot_vel]; 
%                 Session=[Session,data.(u).(s).PSM1.grip_vel]; 
%                 Session=[Session,data.(u).(s).PSM2.pos];
%                 Session=[Session,data.(u).(s).PSM2.rot]; 
%                 Session=[Session,data.(u).(s).PSM2.vel]; 
%                 Session=[Session,data.(u).(s).PSM2.rot_vel]; 
%                 Session=[Session,data.(u).(s).PSM2.grip_vel]; 
% 
%                 k=size(Session,2);
%                 % Adding the gesture
%                 for m=1:length(raw_data.(u).(s).start)
%                     b=raw_data.(u).(s).start(m,1);
%                     e=raw_data.(u).(s).stop(m,1);
%                     Session(b:e,k+1)=raw_data.(u).(s).gestures(m,1);
%                 end
%                 % Cut of the session  according to its gesture
%                 Session=Session(raw_data.(u).(s).start(1):raw_data.(u).(s).stop(end),:);            
%                 if length(find(Session(:,end)==0))~=0
%                     ind=find(Session(:,end)==0);
%                     Session(ind,:)=[];
%                 end
%                 
%                 Features = [Features;Session];
% 
%             end 
%     
%         end