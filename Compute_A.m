function A=Compute_A(x,s); g=9.8;
a42=s.m1*s.ell1*(x(8)*sin(x(2))+x(5)^2*cos(x(2))); a45=2*s.m1*s.ell1*x(5)*sin(x(2));
a43=s.m2*s.ell2*(x(9)*sin(x(3))+x(6)^2*cos(x(3))); a46=2*s.m2*s.ell2*x(6)*sin(x(3));
a52=s.m1*s.ell1*(g*cos(x(2))-x(7)*sin(x(2))); a63=s.m2*s.ell2*(g*cos(x(3))-x(7)*sin(x(3)));
A=[zeros(3) eye(3); 0 -a42 -a43 0 -a45 -a46; 0 a52 0 0 0 0; 0 0 a63 0 0 0];
end % function Compute_A