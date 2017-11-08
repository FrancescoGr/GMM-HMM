 function DATA_SETs = Dataset(User_out,Data_inusers,users)

for i=1:length(User_out)
    u = ['u',int2str(User_out(i))];
    DATA_SETs{i,1}.test_user=users(i); 
    
    for T=1:3
        t = ['t',int2str(T)];
    
        for j=1:numel(fieldnames(Data_inusers.(u)))
             s = ['s',int2str(j)];
             DATA_SETs{i,1}.Test.(t).(s) = Data_inusers.(u).(t).(s);
        end
    end
    
    % Deleting the test user
    train_us=1:length(users);
    train_us(i)=[];
    DATA_SETs{i,1}.Training=[];
    
    
    for y=1:length(train_us)
        u = ['u',int2str(train_us(y))];
        for T=1:3
            t = ['t',int2str(T)];
            if isfield(Data_inusers.(u),(t))
                
                for j=1:numel(fieldnames(Data_inusers.(u).(t)))
                    s = ['s',int2str(j)];
                    DATA_SETs{i,1}.Training = [DATA_SETs{i,1}.Training;Data_inusers.(u).(t).(s)];
                end

                
            end
        
        end

    end
end