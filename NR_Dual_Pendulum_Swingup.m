function [u_k,x_k]=NR_Dual_Pendulum_Swingup(T,u_k)
if nargin<1, T=3; end                 % Default value of T
s.h=0.01; s.N=T/s.h; 
if nargin<2, u_k=zeros(s.N+1,1); end  % Default value of initial u (a pretty bad guess!)
s.mc=10; t=[0:s.N]*s.h; k_max=100;
s.m1=0.2; s.L1=1;    s.ell1=s.L1/2; s.I1=s.m1*s.ell1^2/12;
s.m2=0.1; s.L2=0.5;  s.ell2=s.L2/2; s.I2=s.m2*s.ell2^2/12; alpha=0.1;
s.B=[0; 0; 0; 1; 0; 0]; s.Q=diag([0 0 0 0 0 0]); s.R=0; s.QT=diag([5 40 10 .1 60 10]);
, s.x0=[0; pi; pi; 0; 0; 0]; x_k(1:6,1)=s.x0; res=0;
for k=0:k_max, k
  u=u_k(1); x=s.x0; J=0.25*s.h*(x'*s.Q*x+u'*s.R*u); c=.5; % STEP 1: march/save state
  for n=1:s.N, u=u_k(n);                                  % (from t=0 -> T), compute cost
    f1=RHS(x,u,s); f2=RHS(x+s.h*f1/2,u,s); f3=RHS(x+s.h*f2/2,u,s); f4=RHS(x+s.h*f3,u,s);
    x=x+s.h*(f1/6+(f2+f3)/3+f4/6); x_k(1:6,n+1)=x; u=u_k(n+1);
    x_k(7:9,n)=f1(4:6); if n==s.N, c=.25; end, J=J+c*s.h*(x'*s.Q*x+u'*s.R*u);
  end, f1=RHS(x,u,s); x_k(7:9,s.N+1)=f1(4:6); E=Compute_E(x,s); J=J+0.5*(x'*E'*s.QT*E*x);
  r=s.QT*E*x; g(s.N+1,1)=s.B'*r+s.R*u_k(s.N+1);       % STEPS 2 & 3: march adjoint
  for n=s.N:-1:1, xh=(x_k(:,n+1)+x_k(:,n))/2;         % (from t=T -> 0), compute gradient
    f1=RHSa(r,x_k(:,n+1),s); f2=RHSa(r-s.h*f1/2,xh,s); f3=RHSa(r-s.h*f2/2,xh,s);
    f4=RHSa(r-s.h*f3,x_k(:,n),s); r=r-s.h*(f1/6+(f2+f3)/3+f4/6); g(n,1)=s.B'*r+s.R*u_k(n);
  end, res1=res; res=g'*g;                            % STEPS 4 & 5: update u and repeat
  if (mod(k,4)==0|alpha<1e-4), p_k=-g; else, p_k=-g+p_k*res/res1; end % conjugate gradient
  figure(1); clf; subplot(2,1,1); plot(t,x_k(1,:),'r-',t,x_k(2,:),'b-',t,x_k(3,:),'g-');
                  subplot(2,1,2); plot(t,u_k,'r--');
  [AA,AB,AC,JA,JB,JC]=Bracket(@Compute_J_Ex21_1,0,alpha,J,u_k,p_k,s);  % find triplet
  [alpha,J]=Brent(@Compute_J_Ex21_1,AA,AB,AC,JA,JB,JC,1e-5,u_k,p_k,s)  % refine triplet
  u_k=u_k+alpha*p_k; pause(0.01); if abs(alpha)<1e-12, break, end      % update u_k
end
s.mc=1; for n=1:s.N+1      % Compute u_k corresponding to different s.mc to give same x_k
  E=Compute_E(x_k(1:6,n),s); N=Compute_N(x_k(1:6,n),0,s); u_k(n,1)=s.B'*(E*x_k(4:9,n)-N);
end
end % function NR_Dual_Pendulum_Swingup
