function  Discrepancy=Rentedfinder_HR93_CES(Rented,Pars,n,s)
w=Pars(11);
sigma=Pars(15);
kappa=Pars(16);
gamma=Pars(17);
theta=Pars(5);


    Discrepancy=w*kappa-s*theta*(sigma*n^gamma+(1-sigma)*Rented^gamma)^(theta/gamma-1)*...
        (1-sigma)*Rented^(gamma-1);
    %Discrepancy
end


%MAKE IT WORK



