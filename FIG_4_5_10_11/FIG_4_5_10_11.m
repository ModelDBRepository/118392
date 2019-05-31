%This code generates
%1. Dynamics of a cell-based oscillator (Figure 4, top panel of paper). The
%level(s) of the tonic frive is specified by the user.
%
%2. PRC of a cell-based oscillator if the fast variable is perturbed
%(Figure 4, bottom panel). The adjoint method is used.
%
%3. Coupling function with one-way exitatory connection and mutual
%   inhibition (Figure 5, left and right panel)
%
%4. Coupling function for coupled segments (Figure 10)
%
%5. Stable (decreasing) roots of Figure 10. (Figure 11)
%
disp('Choose values of tonic drive or press return for default [ .05  .07 .09  .11 .13 .15 ]');
Ak=input('Use "[ numbers divided by blank spaces ]" format:' );

if length(Ak)==0
    Ak=[ .05  .07 .09  .11 .13 .15 ]
end



global A th te g szi0 ee;

ee=1;
te=1;
szi0=.1;
g=1.2;



for i=1:length(Ak)
      A=Ak(i);
    disp('Current value of tonic drive:');
    disp(A);
    th=40/(1+(20*A).^2);
    
    %%%%%%%%%dynamics
    [t,y]=period(@wilson_simplified,1000,[.1 0],.00000001,.00000001);
    figure(4);
    hold on;
    subplot(2,length(Ak),i)
    title({'FIG. 4, e=' A});
    plot(2*pi*t/max(t),y); %dynamics of fast and slow variables
    hold on;
    plot(2*pi*t/max(t),NakaRu(y(:,1)-g*y(:,2),szi0),'-'); %mean firing rate of oscillator
    drawnow;

    %%%%%%%%%PRC and coupling function
    [sh{i},H{i}, T(i), fiPRC{i},PRC{i}]=Wilson_sim_coupfn();
    sh{i}=sh{i}/2/pi;

    %plotting PRC
    figure(4);
    subplot(2,length(Ak),length(Ak)+i);
    plot(fiPRC{i},PRC{i}); %PRC
    title({'FIG. 4, e=' A});
end




for i=1:length(H)
    H2{i}=H{i}([(end+1)/2:end 2:(end+1)/2]);%coupling function shifted by \pi
    Hulti{i}=1*H{i}-1*H2{i}; %coupling function of complete segments
    Hinhinh{i}=-H{i}+fliplr(H{i}); %coupling function for mutual inhibition
    
    %plotting coupling function
    figure(5);
    hold on;
    subplot(1,2,1)
    plot(-sh{i}*2*pi,H{i});
    title('FIG 5 left');
    hold on;
    
    %plotting coupling function for mutual inhibition
    subplot(1,2,2)
    title('FIG 5 right');
    plot(2*pi-2*pi*[sh{i}(floor(end/2)+1:end-1) 1+sh{i}(1:floor(end/2)+1)],Hinhinh{i}([floor(end/2)+1:end-1 1:floor(end/2)+1]));
    hold on;
    
    % plotting coupling function for complete segment
    figure(10);
    plot(-sh{i}*2*pi,Hulti{i},'b');
    hold on;
    title('*********FIGURE 10*********');

end

%%%%%%%%%%%%%%%%%%%%%%stable roots of coupling function
Hakt=Hulti; %actual coupling function
for i=1:length(Hakt)
    for j=1:length(Hakt{i})-1
        if Hakt{i}(j+1)*Hakt{i}(j)<=0
            if Hakt{i}(j+1)>Hakt{i}(j)
                stab(i)=(sh{i}(j)*Hakt{i}(j+1)-sh{i}(j+1)*Hakt{i}(j))/(Hakt{i}(j+1)-Hakt{i}(j));
            else
                instab(i)=(sh{i}(j)*Hakt{i}(j+1)-sh{i}(j+1)*Hakt{i}(j))/(Hakt{i}(j+1)-Hakt{i}(j));
            end
        end
    end
end
%plotting roots
figure(11);
plot(1./T(1:end),-stab(1:end),'s-');
title('*********FIGURE 11*********');