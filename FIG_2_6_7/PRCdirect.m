
function [ficycle ycycle fi PRC T]=PRCdirect(funname,x0)

%this function calculates the PRC of the oscillator 'funname' using the direct method. x0
%is an arbitrary initial condition for the integration Input x0 is almost
%arbitrary; 
%
%Parameter numvar in the code determines how many variables areperturbed;
%Parameter tinitial in the code should be chosen such that transients vanish during time t.
%
%Outputs ficycle and ycycle show the dynamics in a cycle;
%Outputs fi, PRC show the phase-response curve
%Output T is the length of the period of the isolated oscillator

tinitial=500;
numvar=3;
%error tolerances
c1=10^(-8);
c2=10^(-8);
%size of perturbation
dy=.001;
%points along the period where the phase rsponse is determined
fi=0:pi/20:2*pi;


% %%%%%%%%ссссссссссссссссссссссссссссссссссссссссссссс
% 
% %%%%%%%%%%% initializing some parameters of @williams
% global ve vl vc exc inh dist ton;
% vl=.01;
% vc=.1;
% % (rows/columns 1,2,...,6 correspond to E,L,C,C,L,E cells , respectively)
% exc0=[0 1 1 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 0]';
% inh0=[0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 1 1 1;1 1 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0]';
% 
% exc=[];
% inh=[];
% ton=[];
% 
% %----- setting some x0 initial condition, vector ton of tonic drives,
% %matrices exc, and inh of excitatory and inhibitory INTRASEGMENTAL connections
% %the effect of intersegmental connections is taken into account in
% %@williams
% x00=[.1 .1 .1 0 0 0]';
% x0=[];
% for i=1:smax
%     x0=[x0;x00];
%     exc=blkdiag(exc,exc0);
%     inh=blkdiag(inh,inh0);
%     ton=[ton; [ve vl vc vc vl ve]'];
% end
% 
% 
% %сссссссссссссссссссссссссссссссссссссссссс

%------------------------length of period and initial point of period
% the oscillator is integrated for some time until transients vanish
[t,y]=ode45(funname,[0 tinitial],x0);
%the oscillator is integrated until the beginning of a new cycle (as
%defined in @newcycle)
for i=1:2
    y0=y(end,:);
    [t,y]=ode45(funname,[0 tinitial],y0,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
end
% length of period and initial point of a cycle are determined
T=t(end);
y0=y(end,:);


%points on one cylce
ficycle=t*2*pi/T;
ycycle=y(:,1:6);

%%%-----------------calculation of the PRC vial integration of oscillators
%%%from unperturbed (y01) and perturbed (y02) initial points

PRC=[];
tphr=fi*T/2/pi; 

for index =1:numvar
    for i=1:length(fi)
        if tphr(i)==0
            t=0;
            y=y0;
        else
            [t,y]=ode45(funname,[0 tphr(i)],y0,odeset('Reltol',c1,'AbsTol',c2));
        end
        y01=y(end,:); %unperturbed initial point
        y02=y01+[zeros(1,index-1) dy zeros(1,length(x0)-index)]; %perturbed initial point
        tphr1=0;
        tphr2=0;

        for j=1:3 %integration through 2-3 cycles
            [t1,y1]=ode45(funname,[0 10*tinitial],y01,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
            [t2,y2]=ode45(funname,[0 10*tinitial],y02,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
            y01=y1(end,:);
            y02=y2(end,:);
            tphr1=tphr1+t1(end);
            tphr2=tphr2+t2(end);
        end
        PRC(i,index)=(tphr1-tphr2)*2*pi/T/dy; %calculation of phase response from time differences
        percon=2*pi/dy;
        PRC=mod(PRC+percon/2,percon)-percon/2; %correction for the case (tphr1-tphr2) accidentally contains a multiple of T 
    end
end

% figure(3);
% plot(fi,PRC);

