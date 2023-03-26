clear all; close all; clc
%% Pendulum Constants
s.mc=10;
s.m1=1; s.L1=2;    s.ell1=s.L1/2; s.I1=s.m1*s.ell1^2/12;
s.m2=.5; s.L2=.5;  s.ell2=s.L2/2; s.I2=s.m2*s.ell2^2/12;

%% Dynamics of Swing Up
s.B=[0; 0; 0; 1; 0; 0]; 
s.B=[0; 0; 0; 1; 0; 0]; s.Q=2*diag([1 1 10 .1 5 20]); s.R=.001; s.QT=2.1*diag([100 240 160 .01 80 60]); % I have not been able to tune this correctly
G = s.B;
Q = s.Q;
R = s.R;
QT = s.QT;

%% Trajectory
T=3;              
s.h=0.01; s.N=T/s.h; 
u_k=zeros(s.N+1,1);  % Default value of initial u (a pretty bad guess!)
t=[0:s.N]*s.h; k_max=100;
alpha=0.1;
s.x0=[0; pi; pi; 0; 0; 0]; x_k(1:6,1)=s.x0; res=0;
for k=0:k_max
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
  figure(1); clf; subplot(2,1,1); plot(t,x_k(1,:),'r-',t,x_k(2,:),'b-',t,x_k(3,:),'g-',t,x_k(5,:),t,x_k(6,:)); legend('x','th1','th2','dth1','dth2');
                  subplot(2,1,2); plot(t,u_k,'r--');
  [AA,AB,AC,JA,JB,JC]=Bracket(@Compute_J_Ex21_1,0,alpha,J,u_k,p_k,s);  % find triplet
  [alpha,J]=Brent(@Compute_J_Ex21_1,AA,AB,AC,JA,JB,JC,1e-5,u_k,p_k,s);  % refine triplet
  u_k=u_k+alpha*p_k; pause(0.01); if abs(alpha)<1e-12, break, end      % update u_k
end
s.mc=1; for n=1:s.N+1      % Compute u_k corresponding to different s.mc to give same x_k
  E=Compute_E(x_k(1:6,n),s); N=Compute_N(x_k(1:6,n),0,s); u_k(n,1)=s.B'*(E*x_k(4:9,n)-N);
end
disp("J = " + J)
disp("States = [" + x_k(1,end) + "," + x_k(2,end) + "," + x_k(3,end) + "," + x_k(4,end) + "," + x_k(5,end) + "," + x_k(6,end) + "]")
%% Control of Swingup (DRE)
K_plus_1 = numel(u_k)+1;
K = zeros(1,6,K_plus_1-1);
eigenvals = zeros(6,6,K_plus_1-1);
X = zeros(6,6,K_plus_1);

for k = numel(u_k):-1:1

    xk = x_k(1:9,k);
    F(:,:,k) = Compute_A(xk,s);
    F_k = F(:,:,k);
    X(:,:,k) = F_k'*X(:,:,k+1)*F_k - F_k'*X(:,:,k+1)*G*pinv(R+G'*X(:,:,k+1)*G)*G'*X(:,:,k+1)*F_k+Q;
    K(1,:,k) = -pinv(R)*G'*pinv(F_k)'*(X(:,:,k)-Q);

end

%% Estimation for Swing Up (DRE)
P0 = diag([1 1 1 1 1 1]);
P = zeros(6,6,K_plus_1);
P(:,:,1) = P0;
L = zeros(6,1,K_plus_1-1);
H = [1 1 1 0 0 0];

for k = 1:K_plus_1-1

    xk = x_k(1:9,k);
    F_k = F(:,:,k);
    P(:,:,k+1) = F_k*P(:,:,k)*F_k' - F_k*P(:,:,k)*H'*pinv(R+H*P(:,:,k)*H')*H*P(:,:,k)*F_k'+Q;
    L(:,1,k) = -P(:,:,k)*H'*pinv(H*P(:,:,k)*H'+R);
    
end
%% Inverted Infinite Horizon Dynamics
syms xdd th1dd th2dd th1 th2 u
g = 9.81;
mc = s.mc; 
m1 = s.m1; l1 = s.L1; I1 = s.I1;
m2 = s.m2; l2 = s.L2; I2 = s.I2;
eqns = [xdd == (u + m1*l1*th1dd + m2*l2*th2dd)*(1/(mc+m1+m2)),th1dd == (m1*g*l1*th1+m1*l1*xdd)*(1/(I1+m1*(l1^2))), th2dd  == (m2*l2*xdd+m2*g*l2*th2)*(1/(I2+m2*(l2^2)))];
solved = solve(eqns,[xdd th1dd th2dd]);
xdd_str = char(solved.xdd);
th1dd_str = char(solved.th1dd);
th2dd_str = char(solved.th2dd);

th1xdd = eqntoss(xdd_str,'th1');
th1th1dd = eqntoss(th1dd_str,'th1');
th1th2dd = eqntoss(th1dd_str,'th1');

th2xdd = eqntoss(xdd_str,'th2');
th2th1dd = eqntoss(th1dd_str,'th2');
th2th2dd = eqntoss(th2dd_str,'th2');

uxdd = eqntoss(xdd_str,'u');
uth1dd = eqntoss(th1dd_str,'u');
uth2dd = eqntoss(th2dd_str,'u');

A_inf = [0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 th1xdd th2xdd 0 0 0; 0 th1th1dd th2th1dd 0 0 0; 0 th1th2dd th2th2dd 0 0 0];
B_inf = [0 0 0 uxdd uth1dd uth2dd]';
C_inf = [1 1 1 0 0 0];
D_inf = 0;
Plant = ss(A_inf,B_inf,C_inf,D_inf,-1);
Plant.InputName = 'un';
Plant.OutputName = 'yt';
%% Control Inverted 
Q_inf = eye(6);
R_inf = .01;
[K_inf,S_inf,P_inf] = lqr(A_inf,B_inf,Q_inf,R_inf); % Uses CARE to solve K

%% Estimation Inverted 
Q_inf = 1;
R_inf = .01;
[kalmf,L_inf,P] = kalman(Plant,Q_inf,R_inf,0);  % Uses CARE to solve K


%% Feedback control Swing Up
% Simulation parameters
t_end = 3;
dt = 0.01;
t = 0:dt:t_end;
num_steps = length(t);

% Preallocate arrays for simulation
x = zeros(9, num_steps);
x_est = zeros(6, num_steps);
y_est = zeros(6, num_steps);
y = zeros(6, num_steps);
u = zeros(1, num_steps);

% Define the disturbance, and noise
w = 0;  % Disturbance
v = 0;  % Measurement noise
reference_signal = 1;
% Initial Condition 
x(2,1) = pi; x(3,1) = pi; 

for k = 1:num_steps-1

    % Plant dynamics
    y(:, k) = C_inf * x(1:6, k) + D_inf * u(:, k) + w;
    F_k = F(:,:,k);
    x(1:6, k+1) = F_k * x(1:6, k) + B_inf * u(:, k) + v;
    
    % Estimator dynamics
    y_est(:,k) = C_inf*x_est(1:6,k);
    x_est(:, k+1) = F_k * x_est(1:6, k) + B_inf * u(:,k) - L(:,1,k)'*(y(:,k)-y_est(:,k));
    
    % Controller dynamics
    u(:, k+1) = K(1,:,k) * x_est(1:6,k) + reference_signal;
end
%% Feedback Control Inverted
for k = num_steps:num_steps*5

    % Plant dynamics
    y(:, k) = C_inf * x(1:6, k) + D_inf * u(:, k) + w;
    x(1:6, k+1) = A_inf * x(1:6, k) + B_inf * u(:, k) + v;
    
    % Estimator dynamics
    y_est(:,k) = C_inf*x_est(1:6,k);
    x_est(:, k+1) = A_inf * x_est(1:6, k) + B_inf * u(:,k) - L_inf'*(y(:,k)-y_est(:,k));
    
    % Controller dynamics
    u(:, k+1) = K_inf * x_est(1:6,k) + reference_signal;
end

%% Plot the results

figure (2)
hold on
subplot(2,1,1)
plot(t, x(1,1:num_steps),t, x(2,1:num_steps),t, x(3,1:num_steps),t, ...
    x(4,1:num_steps),t, x(5,1:num_steps),t, x(6,1:num_steps)); % plot state variables [0,T]
title('Feedback Control of States During Swing Up');
xlabel('Time (s)');
grid on
legend('x','th1','th2','dx_dt','dth1_dt','dth2_dt')
ylabel('State Variables');
subplot(2,1,2)
plot(t,u(1:num_steps))
xlabel('Time (s)');
ylabel('Control Effort');


grid on;
hold off

figure (3)
t = T:.01:T+1205*.01-.01;

subplot(2,1,1)

hold on
plot(t, x(1,num_steps:num_steps*5),t, x(2,num_steps:num_steps*5),...
    t, x(3,num_steps:num_steps*5),t,x(4,num_steps:num_steps*5), ...
    t, x(5,num_steps:num_steps*5),t, x(6,num_steps:num_steps*5)); % plot state variables [T, inf]
legend('x','th1','th2','dx_dt','dth1_dt','dth2_dt')
title('Feedback Control of States During Inverted');
xlabel('Time (s)');
ylabel('State Variables');
subplot(2,1,2)
plot(t,u(num_steps:num_steps*5))
xlabel('Time (s)');
ylabel('Control Effort');


grid on;
hold off

%% Local Functions
function num = eqntoss(char_eqn,var)
    group_len = numel(var);
    firstentry = true;
    thisgroup = false;
    for i = 1:numel(char_eqn)
        if and(char_eqn(i) ~= ' ',firstentry)
            firstentry_i = i;
            firstentry = false;
        end
        if var == char_eqn(i:group_len+i-1) 
           thisgroup = true;
           novar = false;
        end
        if and(char_eqn(i) == ' ',thisgroup)
            lastentry_i = i;
            break
        elseif and(char_eqn(i) == ' ',~thisgroup)
            firstentry = true;
           
        elseif and(i == numel(char_eqn),thisgroup)
            lastentry_i = i;
           
        elseif and(i == numel(char_eqn),~thisgroup)
            novar = true;
        end
    end
    if ~novar
        eval([var '= 1;'])
        eval(['num = ' char_eqn(firstentry_i:lastentry_i) ';'])
    else
        num = nan;
    end

end