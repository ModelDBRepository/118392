%this function is integrated over a period to achieve the coupling function
%output = inp * dEdinp * z
%where
%-  inp = synaptic input (presynaptic activity)
%-  dEdinp is the following derivative: d[right side of fast equation of oscillator] / d inp
%-  z = phase response curve

function output = tolofv(fik,shift)

global fidin ydin fiPRC PRC

for i=1:length(fik)
    fi=mod(fik(i),2*pi);
    y=interp1(fidin,ydin,fi);

    inp=interp1(fidin,ydin(:,1),mod(fi+shift,2*pi));
    z=interp1(fiPRC,PRC,fi);

    dx=10^(-8);
    dEdinp=(wilson_simplified(0,y+[dx 0])-wilson_simplified(0,y-[dx 0]))/2/dx;
    dEdinp=dEdinp(1)+1;
    output(i)=inp*dEdinp*z;
end
