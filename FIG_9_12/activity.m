function [y,ist,dir]=activity(t,x)

% This is an auxiliary function, which helps recording the ends of cycles in
% every segment

y=williams(t,x);
y=y(1+6*[0:length(y)/6-1]);
ist=zeros(size(y));
dir=ist-1;
