function [Data_inusers_rep,LD,K,U,R,T,Gestures]=Subdividing(DatawithStates)

%% Initializations
Gestures = 15; % max(DatawithStates(:,end));
LD= size(DatawithStates,1);
K = 3; %Possible states
U = max(DatawithStates(:,1)); % users B D C
R = max(DatawithStates(:,3)); % possible repetitions
T = max(DatawithStates(:,2)); % possible tasks

start =1;
User(2,:)= [start ,0];

% 1) subdivide each user
for i=2:LD
    if DatawithStates(i-1,1) ~= DatawithStates(i,1)
        u=DatawithStates(i-1,1);
        User(u,:)=[start ,i-1];
        start= i;
    end
end

u=DatawithStates(i-1,1);
User(u,:)=[User(u-1,2)+1 ,i];

for u= 1:U
    if User(u,1)~=0
        Data_inusers{u,1}(:,:)=DatawithStates(User(u,1):User(u,2),:);
    end
end

% 2) subdivide each repetition: Data_inusers_rep(u,r*t)-> t=tasks r= repetitions, u =users
for u= 1:U
    k=1;
    if User(u,1)~=0
        for t= 1:T
            for r= 1:R
                [positions1]=find(Data_inusers{u,1}(:,3)==r);
                [positions2]=find(Data_inusers{u,1}(:,2)==t);
                positions = intersect(positions1,positions2);
                Data_inusers_rep{u,k} = Data_inusers{u,1}(positions,:);
                k=k+1;
            end
        end
        k=1;
    end
end


end