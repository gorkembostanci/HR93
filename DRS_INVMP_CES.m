function  Result=DRS_INVMP_CES(s, w, theta, sigma, kappa,gamma)



rented=((1-sigma)*theta*s/(w*kappa))^(1/(1-theta))*...
    (sigma*(sigma*kappa/(1-sigma))^(gamma/(1-gamma))+1-sigma)^...
    ((theta-gamma)/(gamma*(1-theta)));
Result=(sigma*kappa/(1-sigma))^(1/(1-gamma))*rented;

end