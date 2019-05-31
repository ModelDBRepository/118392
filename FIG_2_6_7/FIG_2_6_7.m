% This code computes certain properties of a 'network-based' segmental
% oscillator model of the lamprey spinal chord proposed by [Williams, 1990].
%
%The following plots are generated:
% 1. the dynamics of an oscillator (Figure 2, top of paper).
% 
% 2. The PRC of an oscillator for perturbations of E,L, and C cells (Figure
%    2, bottom of paper). The direct method is used, i.e. the effect of small finite
%    perturbations is simulated numerically.
%
% 3. Coupling function of two unidirectionally coupled segmental network-based oscillators (Figure 6 of paper), 
%    This is calculated as the integral of 
%    [presynaptic activity] x [difference between reversion potential and current postsynaptic potential] x [postsynaptic PRC] 
%    over a period summed for all connections, at given phaseshift between the connected oscillators
%
%4. the stable (decreasing) roots of the coupling functions (Figure 7 of paper) 
%
%All calculations are done for several values of tonic drive to E cells,
%which can be specified determined by user.


%levels of tonic drive to E cells
disp('Insert level(s) of tonic drive to E cells ');
disp('please write number in [...] divided by blank spaces');
disp('or press return for default [.005 .0075 .01 .02 .04 .06 .07]');
ves=input('>>>:');
if length(ves)== 0
    ves=[.005 .0075 .01 .02 .04 .06 .07];
    disp(ves);
end
%ve=[.005  .01  .07];


%%%%%%%%%%%%%%%%%%%%%%%%%determination of coupling functions;
%%%%elements of Fig 2 are plotted during evaluation of @couplingfuncion
for i=1:length(ves)
    disp('Current level of tonic drive to E cells is')
    disp(ves(i)); disp('Please wait!');
    disp('****************************');
    [shift{i}, H{i}, fi{i}, PRC{i}, T(i)]=couplingfunction(ves(i),i,length(ves));
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%plotting Fig 6
figure(6);
for i=1:length(shift)
    plot(-shift{i},H{i});
    title('***********FIGURE 6*************'); 
   hold on;
end

%%%%%%%%%%%%%%%%%stable zeros of coupling function
disp('calculating roots of copupling functions');
for i=1:length(shift)
    for j=1:length(H{i})-1
        if H{i}(j+1)*H{i}(j)<=0
            if H{i}(j+1)>H{i}(j)
                stab(i)=(shift{i}(j)*H{i}(j+1)-shift{i}(j+1)*H{i}(j))/(H{i}(j+1)-H{i}(j))/2/pi;
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%plotting Fig 7
figure(7);
plot(1./T,-stab,'s-');
title('***********FIGURE 7*************');
disp('Ready!');
