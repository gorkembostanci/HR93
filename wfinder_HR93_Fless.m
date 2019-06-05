function  Discrepancy=wfinder_HR93_Fless(w,v,Pars,Value0)
    Pars(11)=w;
    cE=Pars(3);

    Result=VFI_HR93_Fless(Pars,Value0);
    Value=Result{1};
    
    Discrepancy=(cE-v*(Value(:,1)))^2;
    Discrepancy
end


%MAKE IT WORK



