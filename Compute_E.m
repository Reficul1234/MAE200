function E=Compute_E(x,s); I=eye(3); Z=zeros(3);
E=[I Z; Z [s.mc+s.m1+s.m2         -s.m1*s.ell1*cos(x(2)) -s.m2*s.ell2*cos(x(3));
           -s.m1*s.ell1*cos(x(2))  s.I1+s.m1*s.ell1^2             0            ;
           -s.m2*s.ell2*cos(x(3))          0              s.I2+s.m2*s.ell2^2   ]];
end % function Compute_E