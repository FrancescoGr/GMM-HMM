function [Trans_P] = Transprob(Data_inusers_rep,LD,K,U,R,Gestures)
% Counting the transitions between states within the gestures

% initializations
for i=1:Gestures
    Trans{1,i}=zeros(K,K);
end
Sum=zeros(Gestures,1);

for u= 1:U
    for r= 1:R
        if length(Data_inusers_rep{u,r})~= 0
            LD = length(Data_inusers_rep{u,r});
            Trial=Data_inusers_rep{u,r};
            for i=2:LD

                if Trial(i-1,end-4) == Trial(i,end-4)
                    j= Trial(i,end-4);
                    Trans{1,j}(Trial(i-1,end),Trial(i,end)) = Trans{1,j}(Trial(i-1,end),Trial(i,end))+1;   
                    Sum(j)=Sum(j)+1;
                end

            end
        end
    end
end

for i=1:Gestures
    Trans_P{1,i}=Trans{1,i}/Sum(i);
end


end
