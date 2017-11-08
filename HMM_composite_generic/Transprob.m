function [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,T,Gestures)
% Counting the transitions between states within the gestures

% initializations
for i=1:Gestures
    Trans=zeros(Gestures,Gestures);
end
Sum=0;

for u= 1:U
    k=1;
    for t=1:T
        for r= 1:R
            if length(Data_inusers_rep{u,k})~= 0
                LD = size(Data_inusers_rep{u,k},1);
                Trial=Data_inusers_rep{u,k};
                for i=2:LD
                        Trans(Trial(i-1,end-4),Trial(i,end-4)) = Trans(Trial(i-1,end-4),Trial(i,end-4)) +1;   
                        Sum=Sum+1;           
                end
            end
            k=k+1;
        end
    end
end

Trans_P=Trans./Sum;

end
