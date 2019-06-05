function  Results=GE_HR93_Search(Pars,v)
NGridSize=Pars(12);
SGridSize=Pars(13);
%beta=Pars(1);
A=Pars(2);
%cE=Pars(3);
%cF=Pars(4);
theta=Pars(5);
%tau=Pars(6);
rho=Pars(7);
sigsq_eps=Pars(8);
a=Pars(9);
%p=Pars(10);
w=Pars(11);

entryvector=[zeros(1,SGridSize),v,zeros(1,SGridSize*(NGridSize-2))];


%===============1. Wage Determination         =========================================
tic
options = optimset('Tolfun',1e-2,'MaxFunEvals',10000000,'MaxIter',1000000);
    
    Value0= zeros(SGridSize,NGridSize);
    g=@(Variable)wfinder_HR93(Variable,v,Pars,Value0);

    [w,~]=fsolve(g, w, options);

fprintf('wage determination was done in \n')
toc
%==============2. Calculating Value/Poliy Functions       =============================
Pars(11)=w;
Results=VFI_HR93(Pars,Value0);
Npolicy=Results{2};
Value=Results{1};

%===============3. Grid Construction        ===========================================


[Sgrid, Prob]=mytauchen(a,rho,sigsq_eps,SGridSize);
Sgrid=exp(Sgrid);

NgridLB=DRS_INVMP(Sgrid(1),w,theta);
NgridUB=1.2*DRS_INVMP(Sgrid(SGridSize),w,theta);

Ngrid=linspace(NgridLB,NgridUB,NGridSize);Ngrid(1)=0;


for jj=1:(NGridSize*SGridSize)
n_N(jj)=floor((jj-1)/SGridSize)+1;
n_S(jj)=mod(jj-1,SGridSize)+1;
end

%===============4. Transfer Function        ===========================================

T=zeros(NGridSize*SGridSize,NGridSize*SGridSize);


for ii=(SGridSize+1):(NGridSize*SGridSize)
         a_prime=Ngrid(Npolicy(n_S(ii),n_N(ii)));
        if a_prime==Ngrid(n_N(1))
            B=1;
%         elseif Ngrid(n_N(1))<a_prime && a_prime<Ngrid(n_N(1)+1)
%             B=(a_prime-Ngrid(n_N(1)))/(Ngrid(n_N(2)+1)-Ngrid(n_N(1)));
        else
            B=0;
        end
        for jj=1:SGridSize
            T(ii,jj)=Prob(n_S(ii),n_S(jj))*B;
        end
        if a_prime==Ngrid(n_N(NGridSize*SGridSize))
            B=1;
%         elseif  Ngrid(n_N(NGridSize*SGridSize)-1)<a_prime && a_prime< Ngrid(n_N(NGridSize*SGridSize))
%             B=(a_prime-Ngrid(n_N(NGridSize*SGridSize)-1))/...
%               (Ngrid(n_N(NGridSize*SGridSize))-Ngrid(n_N(NGridSize*SGridSize)-1));
        else
            B=0;
        end
        for jj=1:SGridSize
            T(ii,NGridSize*SGridSize+1-jj)=Prob(n_S(ii),n_S(NGridSize*SGridSize+1-jj))*B;
        end
        
    
    for jj=(SGridSize+1):(NGridSize*SGridSize-1)
        if a_prime==Ngrid(n_N(jj))
            B=1;
%         elseif Ngrid(n_N(jj)-1)<a_prime && a_prime<Ngrid(n_N(jj))
%             B=(a_prime-Ngrid(n_N(jj)-1))/(Ngrid(n_N(jj))-Ngrid(n_N(jj)-1));
%         elseif Ngrid(n_N(jj))<a_prime && a_prime<Ngrid(n_N(jj)+1)
%             B=(a_prime-Ngrid(n_N(jj)))/(Ngrid(n_N(jj)+1)-Ngrid(n_N(jj)));
         else
            B=0;
        end
        T(ii,jj)=Prob(n_S(ii),n_S(jj))*B;
    end
end

lambda=zeros(NGridSize*SGridSize,1)+1/(NGridSize*SGridSize);


%===============5. Entry Determination ===========================================
tic
options = optimset('Tolfun',1e-5,'MaxFunEvals',10000000,'MaxIter',1000000);
    
    M0=0;
    h=@(M)Mfinder_HR93(M, entryvector, lambda, T, NGridSize, SGridSize, A, n_N);

    [M,~]=fsolve(h, M0, options);
    lambda=LambdaCalculator_HR93(M, entryvector, lambda, T, NGridSize, SGridSize, n_N);
fprintf('entry determination was done in \n')
toc


lambda_Matrix=reshape(lambda, [SGridSize,NGridSize]);

%===============6. MisAllocation (Given Entry) ===================================
RealizedOutput=0;

for jj=1:SGridSize
    for kk=1:NGridSize
        RealizedOutput=RealizedOutput+lambda_Matrix(jj,kk)*...
            DRS(Sgrid(jj), Ngrid(kk), theta);
    
    end
end

Results={w,M,lambda,Npolicy, Value, RealizedOutput};
%===============7. Plotting Results ===========================================

figure
for ii=1:SGridSize
subplot (3,SGridSize,ii)
plot(Ngrid,Ngrid(Npolicy(ii,:)),'LineWidth',1);
end
for ii=(SGridSize+1):(2*SGridSize)
subplot (3,SGridSize,ii)
plot(Ngrid,Value(ii-SGridSize,:),'LineWidth',1);
end
subplot (3,SGridSize,2*SGridSize+1)
mesh(Ngrid,Sgrid,lambda_Matrix)

end