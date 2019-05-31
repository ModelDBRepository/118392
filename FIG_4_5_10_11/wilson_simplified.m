
%equation of cell-based oscillator

function y=wilson_simplified(t,x)


global A th te g szi0 ee;

E=x(1);
H=x(2);
y(2)=1/th*(-H+NakaRu(E-g*H,szi0));
y(1)=1/te*(-E+A+ee*NakaRu(E-g*H,szi0));
y=y';



