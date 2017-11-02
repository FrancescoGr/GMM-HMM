function [Groups,States] = Grouping(Data_withGestures,K)
% Grouping of the data, according to the gestures
% K= number of gestures
    for i=1:K
        States{1,i}(:,1)= find( Data_withGestures(:,end) == i);
    end
    for i=1:K
        Groups{1,i}(:,:) = Data_withGestures(States{1,i}(:,1),1:end-1);
    end
end