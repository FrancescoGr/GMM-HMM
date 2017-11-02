function [init] = Startprob(Data_inusers_rep,LD,K,U,R,Gestures)
% Counting of how many initial points for each state

% initializations
for i=1:Gestures
    init(1,i)=0;
end

% Counting which is the first gesture in each trial
for u= 1:U
    for r= 1:R
        if length(Data_inusers_rep{u,r})~= 0
            Gesture_init=Data_inusers_rep{u,r}(1,end-4);
            init(1,Gesture_init) = init(1,Gesture_init)+1;
        end
    end
end

% Definining the probability
Sum=sum(init(1,:));
init=init./Sum;

end