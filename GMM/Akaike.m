function [AIC,BIC] = Akaike(L,p,LT)
% Standard method to evaluate the result of the Algorithm considering the
% high complexity of the problems:

% Akaike criterion
AIC = 2*p+2*(-L);
% AIC_new = 2*p-L;
% Baysian criterion
% BIC = -2*log(L)+p*log(LT);
% BIC_new = L - 0.5*p*log(LT);
% BIC = -2*log(-L)+p*log(LT);
BIC = 2*(-L)+p*log(LT);
end