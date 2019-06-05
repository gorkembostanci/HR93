function  Results=GE_HR93_Hybrid(Pars,v,wupper,wlower)
NGridSize=Pars(12);
SGridSize=Pars(13);
%beta=Pars(1);
A=Pars(2);
cE=Pars(3);
%cF=Pars(4);
theta=Pars(5);
%tau=Pars(6);
rho=Pars(7);
sigsq_eps=Pars(8);
a=Pars(9);
%p=Pars(10);
w=Pars(11);

Value= zeros(SGridSize,NGridSize);


entryvector=[zeros(1,SGridSize),v,zeros(1,SGridSize*(NGridSize-2))];

Discrepancy=1;
Discrepancy1=10;


tic
while abs(Discrepancy)>0.01
    Pars(11)=w;
    Value=VFI_HR93(Pars,Value);
    Value=Value{1};
    Discrepancy=cE-v*(Value(:,1));
    
    if Discrepancy>0
        wupper=w;
        w=(2*wlower+wupper)/3;
    elseif Discrepancy<0
        wlower=w;
        w=(wlower+2*wupper)/3;
    end
    Discrepancy
end
fprintf('wage determination was done in \n')
toc

Results=VFI_HR93(Pars,Value);

Npolicy=Results{2};
Value=Results{1};
[Sgrid, Prob]=mytauchen(a,rho,sigsq_eps,SGridSize);
Sgrid=exp(Sgrid);
%NgridLB=0;
NgridLB=DRS_INVMP(Sgrid(1),w,theta);
NgridUB=DRS_INVMP(Sgrid(SGridSize),w,theta);

Ngrid=linspace(NgridLB,NgridUB,NGridSize);Ngrid(1)=0;
%TensorGrid=kron(Ngrid,Sgrid);

for jj=1:(NGridSize*SGridSize)
n_N(jj)=floor((jj-1)/SGridSize)+1;
n_S(jj)=mod(jj-1,SGridSize)+1;
end

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

%===============ENTRY===========================================
tic
options = optimset('Tolfun',1e-5,'MaxFunEvals',10000000,'MaxIter',1000000);
    
    M0=0;
    h=@(M)Mfinder_HR93(M, entryvector, lambda, T, NGridSize, SGridSize, A, n_N);

    [M,~]=fsolve(h, M0, options);
    lambda=LambdaCalculator_HR93(M, entryvector, lambda, T, NGridSize, SGridSize, n_N);
fprintf('entry determination was done in \n')
toc
%===============================================================
Results={w,M,lambda,Npolicy, Value};


figure
for ii=1:SGridSize
subplot (2,SGridSize,ii)
plot(Ngrid,Ngrid(Npolicy(ii,:)),'LineWidth',1);
end
for ii=(SGridSize+1):(2*SGridSize)
subplot (2,SGridSize,ii)
plot(Ngrid,Value(ii-SGridSize,:),'LineWidth',1);
end

end