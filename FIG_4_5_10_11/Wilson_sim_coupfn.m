% 
%[shift, coupfn, T, fiPRC, PRC, fidin, ydin]=Wilson_sim_coupfn()
%
%This function calculates
%1. the length of period (T) of a cell-based oscillator
%2. the dynamics ydin(fidin) over a period (Figure 4, top panel of paper)
%3. The phase response curve PRC(fiPRC) over a period for perturbations of
%   the fast variable(Figure 4, bottom panel of paper)
%4. The coupling function coupfn(shift) for one-way excitatory coupling
%   between two units (Figure 5, left panel of paper).
%
% The code is optimised for simplicity and not for speed, some computations are evaluated several times. 



function [shift, coupfn, T, fiPRC, PRC, fidin, ydin]=Wilson_sim_coupfn()

global A th te g szi0 ee ec c fidin ydin fiPRC PRC;


dbszam=floor(10*(.05./A).^2);
%The purpose of variable 'dbszam' is the following: the numerical solution of the adjoint problem 
%requires the continuation of an unstable trajectory, the inbuilt ODE solver of matlab is not designed for such purposes.
%To reduce errors, the cycle is cut to dbszam pieces, 
%and the integration is evaluated separately for each piece (from correct initial points).
%The chosen value was found to cause negligible error. Lower values decrease computational need. 
%Too low dbszam is indicated by discontinuous or divergent PRC.

%error tolerance
tolconst=10^(-8);
tolconst2=10^(-8);

%one period
[t,y]=period(@wilson_simplified,1000,[.1 0],tolconst,tolconst2);
T=t(end);

if T==1000
    warning('Not periodic if e=');
    disp(A);
end


fidin=t*2*pi/T;
ydin=y;

y0=y(end,:); %initial point for next integration

%%%%%%%%%%%%%%%%%%% determination of PRC by adjoint method
tkesz=[]; %results are collected here
ykesz=[]; %and here

for i=0:dbszam-1
    if i>0
        [t,y]=ode45(@wilson_simplified,[0 T/dbszam],y0); %integration of oscillator to initial point of i^th piece
        y0=y(end,:);
    end;
    [t1 y1]=ode45(@phaseresponsecurve_Wilson_simplified,[i/dbszam*T T+i/dbszam*T],[1 0 y0]',odeset('Reltol',tolconst,'AbsTol',tolconst2)); %solution of adjoint problem
    [t2 y2]=ode45(@phaseresponsecurve_Wilson_simplified,[i/dbszam*T T+i/dbszam*T],[0 1 y0]',odeset('Reltol',tolconst,'AbsTol',tolconst2)); %solution of adjoint problem

    fit=[y1(size(y1,1),1) y2(size(y2,1),1);y1(size(y1,1),2) y2(size(y2,1),2)]; %fundamental solution matrix at t=T
    [evect evalue]=eig(fit);
    
    %chosing the correct eigenvector
    if abs(evalue(1,1)-1)<abs(evalue(2,2)-1)
        y0phr=(evect*[1 ;0])';
    else
        y0phr=(evect*[0 ;1])';
    end

    %normalisation
    y0phr=2*pi/T*y0phr/(y0phr*wilson_simplified(0,[y0']));
    [t y]=ode45(@phaseresponsecurve_Wilson_simplified,[0 T/dbszam],[y0phr y0]',odeset('Reltol',tolconst));
    
    %adding partial results to tkesz and ykesz
    if i==0 
        tkesz=[tkesz;t+i*T/dbszam];
        ykesz=[ykesz;y];
    else
        tkesz=[tkesz;t(2:end)+i*T/dbszam];
        ykesz=[ykesz;y(2:end,:)];
    end
end

%final normalisation
fiPRC=tkesz/T*2*pi;
PRC=ykesz(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% coupling function
shift=-pi:pi/20:pi; %phase shift between oscillators
for i=1:length(shift)
   coupfn(i) = quad(@(x)tolofv(x,shift(i)),0,2*pi)*T/2/pi; 
end











