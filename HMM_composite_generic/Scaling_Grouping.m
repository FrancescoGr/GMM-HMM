function Data_set=Scaling_Grouping(Data_set,SET)
       
     % Grouping in sets of SET samples (on paper [6], SET = 50),
     % considering every repetition of the task (ring) 

        s=1;
        if SET ~= 1
          % Finding every repetition of the task (ring) 
            first = Data_set(:,end-4);
            second = [first;first(end,end)];
            second = second (2:end,1);

            diff = find(first ~= second);

           % Grouping each repetition in sets of SET
    %         for i=SET:SET:length(Data_set)

            for i=1:length(diff)
                k = diff(i,1);
                if i == 1
                    for j=SET:SET:k
                        Data(s,:) = mean(Data_set(j-SET+1:j,:));
                        s=s+1;
                    end
                else
                    for j=diff(i-1,1)+SET:SET:k
                        Data(s,:) = mean(Data_set(j-SET+1:j,:));
                        s=s+1;
                    end

                end

            end

            clear Data_set
            Data(:,end)=round(Data(:,end));
            Data_set = Data(:,4:end-5);
            Repetition = Data(:,1:3);
        else
            Repetition=Data_set(:,1:3);
            Data = Data_set(:,end-4:end);
            Data_set = Data_set(:,4:end-5);
            
        end
        


      % Standardization
%         for j=1:size(Data_set,2)
%             mu(j)=(1/(size(Data_set,1)))*sum(Data_set(:,j)); % mean of the j-feature
%             standd(j)=std(Data_set(:,j));% std of the j-feature
%         end
% 
%         for i=1:size(Data_set,1)
%             for j=1:size(Data_set,2)
%                 Data_set(i,j)=(Data_set(i,j)-mu(j))/standd(j); % scaling and mean normalization
%             end
%         end


%       % Mean normalization
%         mu=mean(Data_set);
%         MAX=max(Data_set);
%         MIN=min(Data_set);
%         for i=1:size(Data_set,1)
%             for j=1:size(Data_set,2)
%                 Data_set(i,j)=(Data_set(i,j)-mu(1,j))/ (MAX(1,j)-MIN(1,j));
%             end
%         end
%  
        
        Data_set=[Repetition,Data_set,Data(:,end-4:end)];
end