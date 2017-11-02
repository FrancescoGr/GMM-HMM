function [OK,Priors_tot]=StabilityControl(Priors_tot,ComponentProportion,Iter) 

n=size(Priors_tot,1);

if Iter < 100 % Training of the Gaussians to check
    if n < 15 % Creation of the first set of data
        Priors_tot(n+1,:) = sort(ComponentProportion);
    else
        S=std(Priors_tot);       % Gaussians over the data
        M=mean(Priors_tot);
        ComponentProportion=sort(ComponentProportion);
        if (M(1,1)-S(1,1)< ComponentProportion(1,1) &&  ComponentProportion(1,1)< S(1,1)+M(1,1)) && (M(1,2)-S(1,2)< ComponentProportion(1,2) &&  ComponentProportion(1,2)< S(1,2)+M(1,2))&& (M(1,3)-S(1,3)< ComponentProportion(1,3) &&  ComponentProportion(1,3)< S(1,3)+M(1,3))
            Priors_tot(end+1,:) = ComponentProportion;   % since the point belongs to the Gaussians, it can contribute to the dataset
           
            % new_Gaussians over the data
            S_new=std(Priors_tot);
            M_new=mean(Priors_tot);
            
            % Evaluate which point does not belong to the group anymore (cartesian)
            distances = sqrt(sum((Priors_tot-M_new).^2));
            [~,ind]=max(distances);
            
            % Deleting the "ind" data
            Priors_tot(ind,:)=[];
            
        end
    end   
    OK=0;
else % check the new data
    S=std(Priors_tot);       % Gaussians over the data
    M=mean(Priors_tot);
    ComponentProportion=sort(ComponentProportion);
    if (M(1,1)-S(1,1)< ComponentProportion(1,1) &&  ComponentProportion(1,1)< S(1,1)+M(1,1)) && (M(1,2)-S(1,2)< ComponentProportion(1,2) &&  ComponentProportion(1,2)< S(1,2)+M(1,2))&& (M(1,3)-S(1,3)< ComponentProportion(1,3) &&  ComponentProportion(1,3)< S(1,3)+M(1,3))
        OK=1;
    else
        OK=0;
    end
end
end
