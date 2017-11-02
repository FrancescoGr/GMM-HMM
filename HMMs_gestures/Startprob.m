function [init] = Startprob(Data_inusers_rep,LD,K,U,R,Gestures)
% Counting of how many initial points for each state

% initializations
for i=1:Gestures
    init{1,i}=zeros(K,1);
end


% Counting: each time that in a trial, a certain group samples of the same
% gesture starts, we save the particular state.
for u= 1:U
    for r= 1:R
        if length(Data_inusers_rep{u,r})~= 0
            LD = length(Data_inusers_rep{u,r});
            Gesture=Data_inusers_rep{u,r}(1,end-4);
            init{1,Gesture}(Data_inusers_rep{u,r}(1,end),1) = init{1,Gesture}(Data_inusers_rep{u,r}(1,end),1)+1;

            for i=2:LD
                if Data_inusers_rep{u,r}(i-1,end-4) ~= Data_inusers_rep{u,r}(i,end-4)
                    Gesture=Data_inusers_rep{u,r}(i,end-4);
                    init{1,Gesture}(Data_inusers_rep{u,r}(i,end),1) = init{1,Gesture}(Data_inusers_rep{u,r}(i,end),1)+1;
                end
            end
        end
    end
end

% Definining the probability
for i=1:Gestures
    Sum=sum(init{1,i});
    init{1,i}=init{1,i}./Sum;

end
end