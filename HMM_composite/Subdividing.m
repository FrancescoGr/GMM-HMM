function [Data_inusers_rep,LD,K,U,R,Gestures]=Subdividing(DatawithStates)

%% Initializations
Gestures =max(DatawithStates(:,end-4));
LD= size(DatawithStates,1);
K = max(DatawithStates(:,end)); %Possible states
U = max(DatawithStates(:,1)); % users B D C
R = max(DatawithStates(:,2)); % possible repetitions

start =1;
User(1,:)= [start ,0];

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
    Data_inusers{u,1}(:,:)=DatawithStates(User(u,1):User(u,2),:);
end

% 2) subdivide each repetition: Data_inusers_rep(u,r)-> r= repetitions, u =users
for u= 1:U
    for r= 1:R
        [positions]=find(Data_inusers{u,1}(:,2)==r);
        Data_inusers_rep{u,r} = Data_inusers{u,1}(positions,:);
    end
end


end