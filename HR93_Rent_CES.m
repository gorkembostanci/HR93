%---------------------------------------------------------------------------------
%1_Parameters-------------------------------------------------------------
%---------------------------------------------------------------------------------

%1.1_Technology and Preferences-------------------------------------
beta=0.8;        %discount rate
A=0.001;             %disutility of labor

cE=1000;         %entry cost
cF=250;          %fixed cost of operating
theta=0.64;      %drs
sigma=0.4/0.64;  %hired share
gamma=0.75;       %Elasticity of substitution
kappa=1.5;       %renting cost   
    
tau=0.8;         %adj cost parameter

%log(s_t)=a+rho*log(s_{t-1})+epsilon_t

rho=0.93;               %persistence of idiosyncratic shock
sigsq_eps=0.05;         %variance of idiosyncratic shock
a=0.2;                 %constant of  idiosyncratic shock (chosen to generate mean employment of 61.7)

p=1;                %normalized output price

wzero=0.55;        %normalized wage rate (starting point)
Mzero=0.4;          %measure of SS entrants (starting point)

gridamp=1.5;
%1.2_Grids

NGridSize=200;SGridSize=6;

%1.3_Accounting ---------------------------------------------------------------------
Pars=[beta, A, cE, cF, theta, tau, rho, sigsq_eps, a, p, wzero, NGridSize,...
    SGridSize, Mzero, sigma, kappa, gamma, gridamp];
v=zeros(1,SGridSize);
v(4:5)=[0.3,0.7];    %entrant distribution probabilities
ValuePrime_Old= zeros(SGridSize,NGridSize);

mu_zero= zeros(NGridSize,SGridSize)+1/(NGridSize*SGridSize);

%Eqbm - with search================================================================
[Sgrid, Prob]=mytauchen(a,rho,sigsq_eps,SGridSize);
Sgrid=exp(Sgrid);


%-------------Adj Cost----------------------------
tic
Results_Search=GE_HR93_Search_CES(Pars,v);
fprintf('Baseline model solved in \n')
toc

w_GE_Search=Results_Search{1};
M_GE_Search=Results_Search{2};
lambda_GE_Search=reshape(Results_Search{3}, [SGridSize,NGridSize]);
NPolicy_GE_Search=Results_Search{4};
Value_GE_Search=Results_Search{5};  
Output_GE_Search=Results_Search{6};
Rented_GE_Search=Results_Search{7};


NgridLB=DRS_INVMP_CES(Sgrid(1),w_GE_Search, theta, sigma, kappa,gamma);
NgridUB=gridamp*DRS_INVMP_CES(Sgrid(SGridSize),w_GE_Search, theta, sigma, kappa,gamma);

Ngrid=linspace(NgridLB,NgridUB,NGridSize);Ngrid(1)=0;
%-------------Frictionless----------------------------

tic
Results_Fless=GE_HR93_Search_Fless_CES(Pars,v);
fprintf('Frictionless model solved in \n')
toc

w_Fless=Results_Fless{1};
M_Fless=Results_Fless{2};
lambda_Fless=reshape(Results_Fless{3}, [SGridSize,NGridSize]);
NPolicy_Fless=Results_Fless{4};
Value_Fless=Results_Fless{5};  
Output_Fless=Results_Fless{6};
Rented_Fless=Results_Fless{7};


NgridLB=DRS_INVMP_CES(Sgrid(1),w_Fless, theta, sigma, kappa, gamma);
NgridUB=gridamp*DRS_INVMP_CES(Sgrid(SGridSize),w_Fless, theta, sigma, kappa, gamma);

Ngrid_Fless=linspace(NgridLB,NgridUB,NGridSize);Ngrid_Fless(1)=0;
%================================================================================
figure
mesh(Value_GE_Search,'FaceColor','g')
hold on;
mesh(Value_Fless,'FaceColor','r')
hold off;

figure
mesh(Ngrid(NPolicy_GE_Search),'FaceColor','g')
hold on;
mesh(Ngrid_Fless(NPolicy_Fless),'FaceColor','r')
hold off;

figure
mesh(lambda_GE_Search,'FaceColor','g')
hold on;
mesh(lambda_Fless,'FaceColor','r')
hold off;

figure
mesh(Ngrid,Sgrid,Rented_GE_Search,'FaceColor','g'); hold on;
mesh(Ngrid,Sgrid,Rented_Fless,'FaceColor','r'); hold off;
