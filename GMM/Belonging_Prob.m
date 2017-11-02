function [Pb,Pp] = Belonging_Prob(Priors,Covariances,Means,Training_set)
LT = size( Training_set,1);
n = size(Means,2);
     for j = 1:length(Priors)
         for i = 1:LT
             Pb(i,j) = (1/((sqrt((2*pi)^(n)*det(Covariances(:,:,j))))))*exp(-0.5*(Training_set(i,:)-Means(j,:))*pinv(Covariances(:,:,j))*(Training_set(i,:)-Means(j,:))');
         end
     end
     
                 % Posterior Probability (Baysian)
     for i = 1:LT

          Sum = 0;
          for k = 1:length(Priors)
               Sum = Sum+(Pb(i,k)*Priors(1,k));
          end

           for j = 1:length(Priors)
                if Sum == 0
                    l=0;
                end
                Pp(i,j) = (Pb(i,j)*Priors(1,j))/Sum;
           end     
      end  
     
end