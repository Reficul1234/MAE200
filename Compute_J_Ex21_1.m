function J=Compute_J_Ex21_1(u_trial,s);
x=s.x0; u=u_trial(1); J=0.25*s.h*(x'*s.Q*x+u'*s.R*u); c=.5;
for n=1:s.N, u=u_trial(n); if n==s.N, c=.25; end
  f1=RHS(x,u,s); f2=RHS(x+s.h*f1/2,u,s); f3=RHS(x+s.h*f2/2,u,s); f4=RHS(x+s.h*f3,u,s);
  x=x+s.h*(f1/6+(f2+f3)/3+f4/6); J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
end, E=Compute_E(x,s); J=J+0.5*(x'*E'*s.QT*E*x);
end % function Compute_J_Ex21_1