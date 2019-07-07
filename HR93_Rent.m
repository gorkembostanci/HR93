%---------------------------------------------------------------------------------
%1_Parameters-------------------------------------------------------------
%---------------------------------------------------------------------------------

%1.1_Technology and Preferences-------------------------------------
beta=0.8;        %discount rate
A=0.001;             %disutility of labor

cE=1000;         %entry cost
cF=250;          %fixed cost of operating
theta=0.64;      %drs
sigma=0.4;       %hired share
kappa=1.5;       %renting cost   
    
tau=0.8;         %adj cost parameter

%log(s_t)=a+rho*log(s_{t-1})+epsilon_t

rho=0.93;               %persistence of idiosyncratic shock
sigsq_eps=0.05;         %variance of idiosyncratic shock
a=0.2;                 %constant of  idiosyncratic shock (chosen to generate mean employment of 61.7)

p=1;                %normalized output price

wzero=0.55;        %normalized wage rate (starting point)
Mzero=0.1;          %measure of SS entrants (starting point)

%1.2_Grids

NGridSize=200;SGridSize=6;

%1.3_Accounting ---------------------------------------------------------------------
Pars=[beta, A, cE, cF, theta, tau, rho, sigsq_eps, a, p, wzero, NGridSize,...
    SGridSize, Mzero, sigma, kappa];
v=zeros(1,SGridSize);
v(4:5)=[0.3,0.7];    %entrant distribution probabilities
ValuePrime_Old= zeros(SGridSize,NGridSize);

mu_zero= zeros(NGridSize,SGridSize)+1/(NGridSize*SGridSize);


%Eqbm - wage, entry, stats dist================================================


% wupper=1.55;
% wlower=1.51;
% 
% Mupper=0.5;
% Mlower=0.0;
% 
% tic
% Results=GE_HR93(Pars,v,wupper,wlower,Mupper,Mlower);
% toc
% 
% w_GE=Results{1};
% M_GE=Results{2};
% lambda_GE=Results{3};
% NPolicy_GE=Results{4};
% Value_GE=Results{5};  %51 seconds
%================================================================================



%Eqbm - with search================================================================


%-------------Adj Cost----------------------------
tic
Results_Search=GE_HR93_Search(Pars,v);
toc

w_GE_Search=Results_Search{1};
M_GE_Search=Results_Search{2};
lambda_GE_Search=reshape(Results_Search{3}, [SGridSize,NGridSize]);
NPolicy_GE_Search=Results_Search{4};
Value_GE_Search=Results_Search{5};  
Output_GE_Search=Results_Search{6};


%-------------Frictionless----------------------------

tic
Results_Fless=GE_HR93_Search_Fless(Pars,v);
toc

w_Fless=Results_Fless{1};
M_Fless=Results_Fless{2};
lambda_Fless=reshape(Results_Fless{3}, [SGridSize,NGridSize]);
NPolicy_Fless=Results_Fless{4};
Value_Fless=Results_Fless{5};  
Output_Fless=Results_Fless{6};


%================================================================================
figure
mesh(Value_GE_Search,'FaceColor','g')
hold on;
mesh(Value_Fless,'FaceColor','r')
hold off;

figure
mesh(NPolicy_GE_Search,'FaceColor','g')
hold on;
mesh(NPolicy_Fless,'FaceColor','r')
hold off;

figure
mesh(lambda_GE_Search,'FaceColor','g')
hold on;
mesh(lambda_Fless,'FaceColor','r')
hold off;
%Eqbm - Hybrid================================================================


% 
% wupper=1.55;
% wlower=1.51;
% 
% tic
% Results_Hybrid=GE_HR93_Hybrid(Pars,v,wupper,wlower);
% toc
% 
% w_GE_Hybrid=Results_Hybrid{1};
% M_GE_Hybrid=Results_Hybrid{2};
% lambda_GE_Hybrid=Results_Hybrid{3};
% NPolicy_GE_Hybrid=Results_Hybrid{4};
% Value_GE_Hybrid=Results_Hybrid{5}; %51 seconds
%================================================================================

