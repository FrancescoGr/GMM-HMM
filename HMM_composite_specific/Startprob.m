function [init] = Startprob(Data_inusers_rep,LD,K,U,R,T,Gestures)
% Counting of how many initial points for each state

% initializations
for i=1:Gestures
    init(1,i)=0;
end

% Counting which is the first gesture in each trial
for u= 1:U
    k=1;
    for t=1:T
        for r= 1:R
            if length(Data_inusers_rep{u,k})~= 0
                Gesture_init=Data_inusers_rep{u,k}(1,end);
                init(1,Gesture_init) = init(1,Gesture_init)+1;
            end
            k=k+1;
        end
    end
end

% Definining the probability
Sum=sum(init(1,:));
init=init./Sum;

end