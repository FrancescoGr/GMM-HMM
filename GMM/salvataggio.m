clear 
close all
clc


load('Gesture.mat')
load('Covariances.mat')
load('Mean.mat')
load('Pb.mat')
load('Prior.mat')
load('Prob.mat')


for i= 1:length(Gesture)
     gesture=Gesture{1,i};
     if length(gesture)>200
        [Pb{1,i},Pp{1,i}] = Belonging_Prob(Prior{1,i},Covariances{1,i},Mean{1,i},gesture(:,2:end));
        Gesture{1,i}=[Gesture{1,i},Pp{1,i}];
     end
end
save Pp Pp
save Gesture Gesture