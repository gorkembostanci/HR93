function  Discrepancy_M=Mfinder_HR93(M, entryvector, lambda, T, NGridSize, SGridSize, A, n_N, Ngrid)
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
LaborDemand=LaborDemand+lambda(kk)*Ngrid(n_N(kk));
end

Discrepancy_M=abs(LaborDemand-1/A); %SS labor supply equals 1/A in this model
Discrepancy_M
end