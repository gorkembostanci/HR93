function  Result=DRS_INVMP(s, w, theta, sigma, kappa)




Result=((s*sigma)/w)^(1/(1-theta))*((theta-sigma)/(kappa*sigma))^...
    ((theta-sigma)/(1-theta));

end