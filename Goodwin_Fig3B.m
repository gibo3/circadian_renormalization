%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Goodwin model simulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

global k1 k2 k3 p1 p2 r n

%%% Parameters for simulation %%%
tspan = [0 200];        % timespan
dt = 0.001;             % time step size
Enum = 20;              % extract number of frequency

%%% Parameter values %%%
k1 = 0.269;
k2 = 0.200;
k3 = 0.0817;
p1 = 0.290;
p2 = 0.240;
r = 0.180;
n = 15;

%%% Initial values %%%
x0 = [0.756723617646234 0.316943992557932 0.931011055942132];

%%% Simulation for Goodwin model by Runge-Kutta method of order 4 %%%
[t, x] = rungekutta4(@Goodwinmodelf,tspan,x0,dt);

%%% Find peak index %%%
Apeakx = zeros(100000,1);
ipeakx = zeros(100000,1);

dipx1 = 1;
for zi = 2:length(x(:,1))-1
    if x(zi,3) > x(zi-1,3) && x(zi,3) > x(zi+1,3)
        Apeakx(dipx1) = x(zi,3);
        ipeakx(dipx1) = zi;
        dipx1 = dipx1+1;
    end
end

zerox1 = find(Apeakx == 0);
Apeakx(zerox1) = [];
ipeakx(zerox1) = [];

tm1 = ipeakx(1);
tm2 = ipeakx(2);

%%% Period %%%
period = t(tm2)-t(tm1);

%%% Generalized Harmonic Analysis (GHA) %%%
y = x(:,3)-mean(x(tm1:tm2-1,3));
run('GHA.m')

%%% NS value of x3 %%%
NS = sqrt(sum(au)/sum(ad));

%%% phase of second order component, alpha %%%
alpha = mod(-2*sphase(1)+sphase(2),2*pi);
if alpha > pi
    alpha = alpha - 2*pi;
end

%%% Timeseries of x3 %%%
figure
plot(t,x(:,3))
xlim([0 100])
box on

fprintf('period = %4.3f, NS = %4.4f, alpha = %4.4f\n', period, NS, alpha);

%%%%%%%%%%%%%%%%%%%%%
%%% Goodwin model %%%
%%%%%%%%%%%%%%%%%%%%%
function Goodwindef = Goodwinmodelf(~,x)

global k1 k2 k3 p1 p2 r n

x1 = x(1);
x2 = x(2);
x3 = x(3);

x1_dot = r/x3^n-k1*x1;
x2_dot = p1*x1-k2*x2;
x3_dot = p2*x2-k3*x3;

output(1) = x1_dot;
output(2) = x2_dot;
output(3) = x3_dot;

% return a column vector
output = output(:);

Goodwindef = output;
end

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
end