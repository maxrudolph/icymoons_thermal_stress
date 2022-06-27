function beta = ammonia_compressibility(X,P,T)
    delta = .001;
    
    Vc = 1/ammonia_density(X,P,T);
    Vp = 1/ammonia_density(X,P+delta,T);
    Vm = 1/ammonia_density(X,P-delta,T);
    beta = -1/Vc * (Vp-Vm)/(2*delta);

end