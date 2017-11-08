
function [Features,user] = Loader_Gestures()

%% Gestures
     
     user = ['D','C','B','E','G','H','I','F']; 
     sessions= [5,5,5;5,5,5;4,4,5;5,4,5;5,0,5;3,3,4;4,4,5;5,3,5];

    for j=1:length(user)
        u = ['u',int2str(j)];
        for i=1:3
            t = ['t',int2str(i)];
            if i==1
                dir1 = 'G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA_LOUO\JHU-ISI\Training_Set\transcriptions/Knot_Tying_';
                dir2 = 'G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA_LOUO\JHU-ISI\Training_Set/Knot_Tying_';
            end
            if i==2
                dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA_LOUO\JHU-ISI\Training_Set\transcriptions/Needle_Passing_';
                dir2 = 'G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA_LOUO\JHU-ISI\Training_Set/Needle_Passing_';
            end
            if i==3
                dir1='G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA_LOUO\JHU-ISI\Training_Set\transcriptions/Suturing_';
                dir2 = 'G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA_LOUO\JHU-ISI\Training_Set/Suturing_';
            end
            for i_session=1:sessions(j,i)
                if sessions(j,i)~=0
                    s = ['s',int2str(i_session)];

                    % Labels
                    filename = [dir1, user(j),num2str(00),num2str(i_session),'.txt'];
                    raw_data.(u).(t).(s) = ExtractData_gestures(filename);

                    % Kinematic data
                    filename1 = [dir2, user(j),num2str(00),num2str(i_session),'.txt'];
                    data.(u).(t).(s) = ExtractData(filename1);

                    Features.(u).(t).(s)=[];
                    Session=[];
                    Session=[Session,j*(ones(length(data.(u).(t).(s).PSM1.pos),1))];
                    Session=[Session,i*(ones(length(data.(u).(t).(s).PSM1.pos),1))];
                    Session=[Session,i_session*(ones(length(data.(u).(t).(s).PSM1.pos),1))];   
                    Session=[Session,data.(u).(t).(s).PSM1.pos];
                    Session=[Session,data.(u).(t).(s).PSM1.rot]; 
                    Session=[Session,data.(u).(t).(s).PSM1.vel]; 
                    Session=[Session,data.(u).(t).(s).PSM1.rot_vel]; 
                    Session=[Session,data.(u).(t).(s).PSM1.grip_vel]; 
                    Session=[Session,data.(u).(t).(s).PSM2.pos];
                    Session=[Session,data.(u).(t).(s).PSM2.rot]; 
                    Session=[Session,data.(u).(t).(s).PSM2.vel]; 
                    Session=[Session,data.(u).(t).(s).PSM2.rot_vel]; 
                    Session=[Session,data.(u).(t).(s).PSM2.grip_vel]; 

                    k=size(Session,2);

                    % Adding the gesture
                    for m=1:length(raw_data.(u).(t).(s).start)
                        b=raw_data.(u).(t).(s).start(m,1);
                        e=raw_data.(u).(t).(s).stop(m,1);
                        Session(b:e,k+1)=raw_data.(u).(t).(s).gestures(m,1);
                    end

                    % Cut of the session  according to its gesture
                    Session=Session(raw_data.(u).(t).(s).start(1):raw_data.(u).(t).(s).stop(end),:);            
                    if length(find(Session(:,end)==0))~=0
                        ind=find(Session(:,end)==0);
                        Session(ind,:)=[];
                    end

                    Features.(u).(t).(s) = Session;
                end

            end 
    
        end
    end

%      sessions= [1:15];
%      user = ['D','C','B']; 
%      n=1;
% 
%     for j=1:length(user)
%         u = ['u',int2str(j)];
%         for i_session = sessions
%             s = ['s',int2str(i_session)];
%             filename = ['C:\Users\Francesco-Greg\Documents\MATLAB\HMM\DATA\JHU-ISI\Training_Set\transcriptions/',user(j),num2str(00),num2str(i_session),'.txt'];
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
%                 filename = ['C:\Users\Francesco-Greg\Documents\MATLAB\HMM\DATA\JHU-ISI\Training_Set/',user(j),num2str(00),num2str(i_session),'.txt'];
%                 data.(u).(s) = ExtractData(filename);
%             end
%         end
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
%                 Session=[Session,i_session*(ones(length(data.(u).(s).PSM1.pos),1))];                Session=[Session,data.(u).(s).PSM1.pos];
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
end