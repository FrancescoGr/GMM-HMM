function L = Evaluation(LT,Pb,Priors)

        S = log(Pb*Priors);
        L = sum(S);

end
