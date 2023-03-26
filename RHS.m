function R=RHS(x,u,s);  E=Compute_E(x,s); N=Compute_N(x,u,s); R=E\N;
end % function RHS