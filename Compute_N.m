function N=Compute_N(x,u,s);
N=[x(4); x(5); x(6); -s.m1*s.ell1*sin(x(2))*x(5)^2-s.m2*s.ell2*sin(x(3))*x(6)^2+u;
 s.m1*9.8*s.ell1*sin(x(2)); s.m2*9.8*s.ell2*sin(x(3)) ];
end % function Compute_N