function dydt = full_ode_system(t,y,k1,k2,kp3,km3,eT,s0)

s = y(1);
e = y(2);
cs = y(3);
i = y(4);
ci = y(5);

dsdt = -k1*e*eT*s;
dedt = -k1*e*s*s0+k2*cs-kp3*s0*e*i+km3*ci;
dcsdt = k1*e*s*s0-k2*cs;
didt = -kp3*e*i*eT + km3*ci*eT/s0;
dcidt = kp3*e*i*s0-km3*ci;

dydt = [dsdt; dedt; dcsdt; didt; dcidt];
    
end