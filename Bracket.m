function [AA,AB,AC,JA,JB,JC]=Bracket(Compute_J,AA,AB,JA,X,P,V)
% INPUT: {AA,AB} are guesses of A near a minimum of J(A)=Compute_J(X+A*P), with JA=J(AA).
% OUTPUT: {AA,AB,AC} bracket the minimum of J(A), with values {JA,JB,JC}.
JB=Compute_J(X+AB*P,V);  if JB>JA; [AA,AB]=Swap(AA,AB); [JA,JB]=Swap(JA,JB); end
AC=AB+2*(AB-AA); JC=Compute_J(X+AC*P,V);
while (JB>JC)                                 % At this point, {AA,AB,AC} has JA>JB>JC.
  AN=AC+2.0*(AC-AB);  JN=Compute_J(X+AN*P,V);  % Compute new point AN outside of triplet
  AA=AB; AB=AC; AC=AN;  JA=JB; JB=JC; JC=JN;  % {AB,AC,AN} -> {AA,AB,AC}
end
end % function Bracket