%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Runge-Kutta method of order 4 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% f     : function
% x0    : initial value
% tspan : time span
% dt    : discrete time size
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, x] = rungekutta4(f, tspan, x0, dt)
t = (tspan(1):dt:tspan(2))';
x = zeros(size(t,1),size(x0,2));
x(1,:) = x0;

for l = 1:size(t,1)-1
    k1 = f(1,x(l,:))*dt;
    k2 = f(1,x(l,:)+k1'/2)*dt;
    k3 = f(1,x(l,:)+k2'/2)*dt;
    k4 = f(1,x(l,:)+k3')*dt;
    
    x(l+1,:) = x(l,:)+(k1+2*k2+2*k3+k4)'/6;
end
return;