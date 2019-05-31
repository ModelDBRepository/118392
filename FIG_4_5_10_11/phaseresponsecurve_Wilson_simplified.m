function y=phaseresponsecurve_Wilson_simplified(t,x)

%x(1,2) adjoint problem
%x(3,4) original problem 


y=linearizedWilson_simplified(x); %y(1:2) is needed
y(3:4)=wilson_simplified(0,x(3:4)); %this call overwrites y(3:4)