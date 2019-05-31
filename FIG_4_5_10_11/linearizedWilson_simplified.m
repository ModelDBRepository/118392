


function y=linearizedWilson_simplified(x)

global A th te g szi0 ee;

xsajat=x(1:2);
xorb=x(3:4);

y=-DF_simplified(x(3),x(4))'*xsajat; %adjoint problem: y=-Df'*x;

