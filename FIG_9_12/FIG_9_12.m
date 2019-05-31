%This code performs simulation of a chain of smax network-oscillators with local (but not only nearest-neighbor)
%connectivity. Outputs are:
%1. Phase lags along the chain (Figure 12)
%2. Average phase lag along the chain (one datapoint of Figure 9)


%--------these global variables are passed to @williams 
clear;

global ve vl vc exc inh dist ton



%number of segmental oscillators in chain
smax=30;

%--------Tonic drive to E,L,C cells

ve=input('Select level of tonic drive to E cells (between 0.005 and 0.07, default=0.02): ');
if length(ve)==0
    ve=0.02;
end

vl=.01;
vc=.1;
%--------parameters of intersegmental connection strength
%the strength of length d connections is e0*e^d for ascending connections
%and h0*h^d for descending ones (or 0 if d>5)
e=.5;
h=.5;
e0=1;
%h0 comes as user input
disp('Please type the strength of descending connections relative to ascending ones!')
disp('suggested values: 0 for unidirectional coupling, 0.2 for bidirectional')
h0=input('Current value:  ');
%matrix of connection strengths
dist=[-5 -4 -3 -2 -1 0 1 2 3 4 5;h0*[h^5 h^4 h^3 h^2 h^1] 1 e0*[e e^2 e^3 e^4 e^5]]';

%--------------------------------
% connectivity matrix of exitatory and inhibitory connections in a segments
% according to Fig. 1 (without E-E connections)


% (rows/columns 1,2,...,6 correspond to E,L,C,C,L,E cells , respectively)
exc0=[0 1 1 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 0 0 0;0 0 0 1 1 0]';
inh0=[0 0 0 0 0 0;0 0 1 0 0 0;0 0 0 1 1 1;1 1 1 0 0 0;0 0 0 1 0 0;0 0 0 0 0 0]';

exc=[];
inh=[];
ton=[];

%----- setting some x0 initial condition, vector ton of tonic drives,
%matrices exc, and inh of excitatory and inhibitory INTRASEGMENTAL connections
%the effect of intersegmental connections is taken into account in
%@williams
x00=[.1 .1 .1 0 0 0]';
x0=[];
for i=1:smax
    x0=[x0;x00];
    exc=blkdiag(exc,exc0);
    inh=blkdiag(inh,inh0);
    ton=[ton; [ve vl vc vc vl ve]'];
end

%-----------integration of ODE for some time and recording ends of cycles
% in variable e; End of cycle is found when the left E cell in a segment
% crosses 0 
% e is a list of times when any of the smax oscillators begins a new cycle
[t,y,e]=ode45(@williamschain, [0,800],x0,odeset('events',@activity));

%calculation of period from e
T=(sum(e(end-smax+1:end))-sum(e(end-2*smax+1:end-smax)))/smax;

%time-lags in e normalised by T
dif=diff(e)/T;

%calculating mean value m of phase lags between neighbors in last 25% of
%time (to filter out transients)
m=[];
for j=1:smax
m(j)=mean(dif(j+smax*[floor(length(dif)/smax*3/4):-1+floor(length(dif)/smax)]));
end

%cyclic permutation of m to begin with lag between segments 1 and 2
[a ii]=max(m);
m=m([ii+1:end 1:ii]);

%%%%%%%%------------------plotting figure
figure(12);
title('**********FIGURE 12***********');
hold on;
plot(m(1:end-1),'o-'); 

figure(9);
hold on;
title('Data points in Figure 9')
plot(1/T,mean(m(1:end-1)),'o');

% figure();
% plot(t(200:500),y(200:500,1:3));