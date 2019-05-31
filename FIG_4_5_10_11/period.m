%integration of one period of cell-based oscillator
function [t,y]=period(fuggveny,tmax,y0,c1,c2)

for i=1:4
    [t,y]=ode45(fuggveny,[0 tmax],y0,odeset('events',@newcycle,'Reltol',c1,'AbsTol',c2));
    y0=y(end,:);
end

