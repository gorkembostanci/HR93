function  Result=AdjCostHR(tau, n, nprime)

%Result=tau*max(0,n-nprime);
Result=tau*abs(n-nprime);
%Result=tau*(abs(n-nprime))^1.1;
%Result=tau*abs(n-nprime)+200*tau*any(abs(n-nprime));
end