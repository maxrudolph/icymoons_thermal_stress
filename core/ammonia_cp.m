function Cp = ammonia_cp(X)
    Cp_h20 = 4182; %J/kg/K
    Cp_Nh3 = 4744; %J/kg/K
    Cp = X*Cp_Nh3 + (1-X)*Cp_h20;
end