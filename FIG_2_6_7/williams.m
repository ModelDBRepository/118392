function y=williams(t,x)

ve=x(7);
vl=x(8);
vc=x(9);
%parameter=x(7);
x=x(1:6);



VIN=[0 1 1 0 0 0;0 0 -1 0 0 0;0 0 0 -1 -1 -1];
VIN=[VIN;-1 -1 -1 0 0 0;0 0 0 -1 0 0;0 0 0 1 1 0];


%exc=excszorzo*[0.025 0.01 .1 .1 0.01 0.025]'; %excitation from brain
exc=[ve vl vc vc vl ve]';
%E-L-C-C-L-E
th=0; %activity threshold
tau=.1;%timescaling
tcells=[1 1 1 1 1 1]';
r=0;%resting potential
vex=1;
W=abs(sign(VIN));

% switch parameter
%     case {1}


X=x*ones(1,length(x)); %x of other unit; X' used for x ofown unit
% VIN([2 5],:)=0;
% VIN(:,[2 5])=0;
 %connection weights all 1
%exc=exc+x-th;%overall activity of cells
% W
% (X-th)
% (VIN-X')
% % (sum(W.*(X-th).*(VIN-X'),1))'
y=tcells.*(exc.*(vex-x)+tau*(r-x)+(sum(W.*max(0,X-th).*(VIN-X'),1))');

y(7:9)=0;


%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %feloszcillator
% X=x*ones(1,length(x)); %x of other unit; X' used for x of own unit
% exc=excszorzo*[.1 .01 .1 ]'; %excitation from brain
% %E-L-C
% th=0; %activity threshold
% tau=.1;%timescaling
% r=0;%resting potential
% vex=1;
% tcells=[1 1 1]';
% VIN=[0 1 1 ;0 0 -1 ;0 0 0 ];
% W=abs(sign(VIN)); %connection weights all 1
% %exc=exc+x-th;%overall activity of cells
% % W
% % (X-th)
% % (VIN-X')
% % (sum(W.*(X-th).*(VIN-X'),1))'
% y=tcells.*(exc.*(vex-x)+tau*(r-x)+(sum(W.*(max(0,X-th)).*(VIN-X'),1))');