function  Result=DRS_CES(s, n, theta, rented, sigma, gamma)

Result=s*(sigma*n^gamma+(1-sigma)*rented^gamma)^(theta/gamma);


end