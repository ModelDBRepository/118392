function J=DF_simplified(E,H) 
%JAcobian of the cell-based CPG model eq.-s

global A th te g szi0 ee;


dx=10^(-8);
J(:,1)=(wilson_simplified(0,[E+dx H])-wilson_simplified(0,[E-dx H]))/2/dx;
J(:,2)=(wilson_simplified(0,[E H+dx])-wilson_simplified(0,[E H-dx]))/2/dx;



