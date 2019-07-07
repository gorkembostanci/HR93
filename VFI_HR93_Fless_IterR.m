function  Results=VFI_HR93_IterR(Pars,ValuePrime_Old)
NGridSize=Pars(12);
SGridSize=Pars(13);
beta=Pars(1);
%cE=Pars(3);
cF=Pars(4);
theta=Pars(5);
tau=Pars(6);
rho=Pars(7);
sigsq_eps=Pars(8);
a=Pars(9);
p=Pars(10);
w=Pars(11);
sigma=Pars(15);
kappa=Pars(16);


%--------------------2-Initialize -------------------------

ValuePrime=ValuePrime_Old;
Npolicy= zeros(SGridSize,NGridSize);
%--------------------3-Deterministic Steady State and Grid Formation --
[Sgrid, Prob]=mytauchen(a,rho,sigsq_eps,SGridSize);
Sgrid=exp(Sgrid);
%Sgrid(1)=0;
%NgridLB=0;
NgridLB=DRS_INVMP(Sgrid(1),w,theta, sigma, kappa);
NgridUB=DRS_INVMP(Sgrid(SGridSize),w,theta, sigma, kappa);

Ngrid=linspace(NgridLB,NgridUB,NGridSize);
%Ngrid2=logspace(log(NgridLB),log(NgridUB),NGridSize);
Ngrid(1)=0;
%Ngrid2(1)=0;
%Ngrid(250)
%--------------------4-Iteration -------------------------
%options = optimset('Tolfun',1e-8,'MaxFunEvals',10000000,'MaxIter',1000000);


Discrepancy=1;
while Discrepancy>0.000001
    for Sstate = 1:SGridSize   %Iterating over states
        Nchosen=1;
        for Nstate = 1:NGridSize
            %Utilize monotonicity and concavity later on
            ContUtilPrime=-AdjCostHR_Fless(tau, Ngrid(Nstate), 0);
            star=0;
            for Nchoice=Nchosen:NGridSize   %Iterating over choices
                
                RentedChoice=(((theta-sigma)*Sgrid(Sstate)*Ngrid(Nchoice)^sigma)/...
                    (w*kappa))^(1/(1-theta+sigma));
                
                Profit= p*DRS(Sgrid(Sstate), Ngrid(Nchoice), theta, RentedChoice, sigma)-...
                    w*Ngrid(Nchoice)-w*kappa*RentedChoice-...
                    AdjCostHR_Fless(tau, Ngrid(Nstate), Ngrid(Nchoice));
                FutureUtil=0;
                
                for Snext=1:SGridSize
                    FutureUtil = FutureUtil+...
                        Prob(Sstate,Snext)*ValuePrime_Old(Snext,Nchoice);
                end                
                                        
                ContUtil=Profit+beta*FutureUtil-p*cF;
                
                if Nchoice==1
                    ContUtil=-AdjCostHR_Fless(tau, Ngrid(Nstate), 0);
                end
                
                if ContUtil>ContUtilPrime
                    ContUtilPrime=ContUtil;
                    Nchosen=Nchoice;
                    star=1;
                elseif star==1 && ContUtil<ContUtilPrime
                    break
%                     fprintf('broken at ')
%                     Nchoice
                end                
            end
            ValuePrime(Sstate,Nstate)=ContUtilPrime;
            Npolicy(Sstate,Nstate)=Nchosen;
            
        end        
    end
    
    Discrepancy=max(abs(ValuePrime-ValuePrime_Old));
    ValuePrime_Old=ValuePrime; 
    Discrepancy2=1;
    
    while Discrepancy2>0.000001
        for Sstate = 1:SGridSize   %Iterating over states
            for Nstate = 1:NGridSize
                %Utilize monotonicity and concavity later on
                %star=0;
                RentedChoice=(((theta-sigma)*Sgrid(Sstate)*Ngrid(Npolicy(Sstate,Nstate))^sigma)/...
                    (w*kappa))^(1/(1-theta+sigma));
                
                Profit= p*DRS(Sgrid(Sstate), Ngrid(Npolicy(Sstate,Nstate)), theta, RentedChoice, sigma)-...
                    w*Ngrid(Npolicy(Sstate,Nstate))-w*kappa*RentedChoice-...
                    AdjCostHR_Fless(tau, Ngrid(Nstate), Ngrid(Npolicy(Sstate,Nstate)));
                FutureUtil=0;
                
                for Snext=1:SGridSize
                    FutureUtil = FutureUtil+...
                        Prob(Sstate,Snext)*ValuePrime_Old(Snext,Npolicy(Sstate,Nstate));
                end 
                
                ContUtilV=Profit+beta*FutureUtil-p*cF;
                
                if Npolicy(Sstate,Nstate)==1
                    ContUtilV=-AdjCostHR_Fless(tau, Ngrid(Nstate), 0);
                end
                ValuePrime(Sstate,Nstate)=ContUtilV;
            end
            
        end
        
        
        
        Discrepancy2=max(abs(ValuePrime-ValuePrime_Old));
        ValuePrime_Old=ValuePrime;
    end
    
end

Results={ValuePrime, Npolicy};

% figure
% for ii=1:6
% subplot (4,3,ii)
% plot(Ngrid,Ngrid(Npolicy(ii,:)),'LineWidth',1);
% end
% for ii=7:12
% subplot (4,3,ii)
% plot(Ngrid,ValuePrime(ii-6,:),'LineWidth',1);
% end

end