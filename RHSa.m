function R=RHSa(r,x,s); E=Compute_E(x,s); A=Compute_A(x,s); R=-E'\(A'*r+s.Q*x(1:6));
end % function RHSa