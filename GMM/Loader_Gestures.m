
function Features = Loader_Gestures()

%% Gestures
     
     user = ['D','C','B','E','G','H','I']; 
     sessions= [15;15;13;14;10;10;13];

    for j=1:length(user)
        u = ['u',int2str(j)];
        for i_session =1: sessions(j)
            s = ['s',int2str(i_session)];
            filename = ['G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Training_Set\transcriptions/',user(j),num2str(00),num2str(i_session),'.txt'];
            raw_data.(u).(s) = ExtractData_gestures(filename);
        end
    end


    %% Pairing Gestures to real data
        for j=1:length(user)
        u = ['u',int2str(j)];
            for i_session =1: sessions(j)
                s = ['s',int2str(i_session)];
                filename = ['G:\Politecnico\V annob\Thesis\HMM\CODICE\DATA\JHU-ISI\Training_Set/',user(j),num2str(00),num2str(i_session),'.txt'];
                data.(u).(s) = ExtractData(filename);
            end
        end

    % Features extraction (useful attributes)
    Features=[];

        for j=1:length(user)
            u = ['u',int2str(j)];
            for i_session  =1: sessions(j)

                s = ['s',int2str(i_session)];
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
                
                Features = [Features;Session];

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