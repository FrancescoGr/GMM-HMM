function [errorM,errorC,errorP] = Err_conv(Means_test,check,Covariances_test,Priors_test)

if check ~=1
    for i = 1: length(Means_test)-1
            errorMt(i,1)=sum(rms(Means_test{1,check}-Means_test{1,i}));
            errorCt(i,1)=sum(sum(rms(Covariances_test{1,check}-Covariances_test{1,i})));
            errorPr(i,1)=rms(Priors_test{1,check}-Priors_test{1,i});
    end
     
    [errorM,indM] = min(errorMt);
    [errorC,indC] = min(errorCt);
    [errorP,indP] = min(errorPr);
    
    if indM ~= indC || indC ~= indP 
        errorM=10;
        errorC=10;
        errorP=10;
    end
else
    errorM=10;
    errorC=10;
    errorP=10;

end
    
end