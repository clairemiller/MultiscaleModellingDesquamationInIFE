function [k1,k2] = cnd_params(pH, klk)
    if (klk == 7)
        k2 = 0.601*10^3;
        k1 = 7.4945 * 10^7;
    elseif (klk == 5)
        k2 = 2.2849*10^3;
        k1 = 4.9703*10^7;
    else
        error("Only value for klk inputs of 5 or 7")
    end
    
end