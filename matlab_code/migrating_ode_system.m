function dydt = migrating_ode_system(t,y,v,eT,iT,s0,klk)
    z = v*t;
    pH = 6.8482 - 0.3765*z - 5.1663*z^2 + 3.1792*z^3;
    [k1,k2] = cnd_params(pH,klk);
    [kp3,km3] = lekti_params(pH, klk);

    dydt = full_ode_system(t,y,k1,k2,kp3,km3,eT,s0);
end