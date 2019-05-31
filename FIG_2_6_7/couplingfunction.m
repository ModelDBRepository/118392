%[shift H fi PRC T]=couplingfunction(ve,i,j)
%
%this function determines: 
%
%-the length of the period T of a network based oscillator; 
%-the phase response curves PRC(fi) 
%-the coupling function H(shift)
%
%if tonic drive to E cells is ve; i,j determine the index of the subplot
%where the results are displayed and the total number of subplots.
%
%the dynamics of an oscillator (Fig 2) is also plotted 


function [shift H fi PRC T]=couplingfunction(ve,i,j)



%%%%%%%%%%%%%%%%%%%%% PRC
[fi0 y0 fi PRC T]=PRCdirect(@williams,[.1 0 0 0 0 0 [ve .01 .1]]);
y0=y0(:,1:6);
fi=fi';
%%%%%%%%%%%%%%%%%%%%% plotting the dynamics of an oscillator
figure(2);
if i==1
    title ('************FIGURE 2*************');
end
subplot(2,j,i);

hold on;
plot(fi0,y0(:,1:3));
hold off;
subplot(2,j,j+i);
plot(fi,PRC);
title({'FIG. 2, tonic drive =' ve});
hold off;


%%%%%%%%%%%%%%%%%%%%
dfi=fi(2)-fi(1);
dt=dfi/(2*pi/T);
y=interp1(fi0,y0,fi,'spline');
shift=-pi:pi/20:pi;%phase-shifts at which the coupling function is determined

%connectivity matrix of segmental oscillators. Rows 1,2,...,6 mean
%E,L,C,C,L,E cells, repectively; +1=excitatory connection; -1=inhibitory
%connection
conn=[0 0 0 -1 0 0; 1 0 0 -1 0 0;1 -1 0 -1 0 0];  
conn=[conn;0 0 0 0 0 0;1 0 0 0 0 0 ; 0 0 0 -1 0 0];
conn=[conn;flipud(fliplr(conn))];

% numerical integration of coupling function
pert=[];
H=[];
for k=1:length(shift) %phase shift between coupled oscillators
    for i=1:3 %the effect of incoming connections of cell i
        for j=1:length(fi)-1 %integrating over a period
            w=abs(conn(i,:));%outgoing connection strengths of i^th cell 
            
            fr=interp1(fi',y,mod(fi(j)+shift(k),2*pi),'spline'); %calculating firing rates
            fr=max(fr,0);
            
            va=sign(conn(i,:))-y(j,i); %this term takes inzto account that connections have reversion potential (+1 or -1)
            pert(j,i)=sum(fr.*w.*va*dt); %perturbations from incoming connections
        end
    end

    %The coupling function is determined from the integral of pert*PRC over
    %a period summed for all elements of the oscillator.
    h=pert.*PRC(1:end-1,:);
    H(k)=2*sum(sum(h)); %'2*' because only perturbations on half of the segment were integrated
end






