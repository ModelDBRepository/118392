function y=williams(t,x)
global ve vl vc exc inh dist ton;

%number of segmental units in chain
smax=length(x)/6;

%some 0-s are added to vector x for technical reasons
xx=[zeros(6*5,1); x; zeros(6*5,1)];



%dynamics of uncoupled oscillator
y=ton.*(1-x)-.1*x;

for i=1:size(dist,1)
    y=y+dist(i,2)*exc*max(0,xx(31+dist(i,1)*6:end-30+dist(i,1)*6)).*(1-x);
    y=y+dist(i,2)*inh*max(0,xx(31+dist(i,1)*6:end-30+dist(i,1)*6)).*(-1-x);
end
