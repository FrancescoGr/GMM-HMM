function Beginning_P = Startprob(DatawithStates,LD,K);
% Counting of how many initial points for each state
init=zeros(K,1);
init(DatawithStates(1,end),1)=1;

for i=2:LD
    if DatawithStates(i-1,1) ~= DatawithStates(i,1)
        for j=1:K
           if(DatawithStates(i,end) == j)
               init(j,1) = init(j,1)+1;
           end
        end
    end
end

Beginning_P = init./sum(init);

end