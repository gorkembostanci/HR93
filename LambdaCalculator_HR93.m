function  Results_lambda=LambdaCalculator_HR93(M, entryvector, lambda, T, ...
    NGridSize, SGridSize, n_N)
entry=M*entryvector;
Discrepancy2=1;

while abs(Discrepancy2)>0.0001
    
   lambdaprime=T'*lambda+entry';
   Discrepancy2=max(abs(lambda-lambdaprime));
   lambda=lambdaprime;
   %Discrepancy2
end
LaborDemand=0;
for kk=1:NGridSize*SGridSize
LaborDemand=LaborDemand+lambda(kk)*n_N(kk);
end

Results_lambda=lambda;

end