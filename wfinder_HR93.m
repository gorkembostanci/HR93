function  Discrepancy=wfinder_HR93(w,v,Pars,Value0)
    Pars(11)=w;
    cE=Pars(3);

    Result=VFI_HR93_IterR(Pars,Value0);
    Value=Result{1};
    
    Discrepancy=(cE-v*(Value(:,1)))^2;
    %Discrepancy
end


%MAKE IT WORK



