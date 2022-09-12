function [kp3, km3] = lekti_params(pH, klk)
    if (klk == 7)
        kp3 = 7e3*exp(1.8*pH);
        km3 = 2.9e4*exp(-1.1*pH);
    elseif (klk == 5)
        kp3 = (5.2*pH-19.5)*10^7;
        km3 = 2.3e6*exp(-3*pH);
    else
        kp3 = 0.5*( 7e3*exp(1.8*pH) + (5.2*pH-19.5)*10^7 );
        km3 = 0.5*( 2.9e4*exp(-1.1*pH) + 2.3e6*exp(-3*pH) );
    end    
end