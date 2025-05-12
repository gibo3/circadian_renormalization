%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 4 variables model calculation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

datetime('now')

global k1 k2 k3 p1 p2 v K n L Omega

%%% Parameters for simulation %%%
tspan = [0 200];        % timespan
dt = 0.001;             % time step size

%%% Activation energy and frequency factor %%%
% Activation energy
Ei = zeros(1,6);
Ei(1) = 1.39392e+04;    % k1
Ei(2) = 6.31343e+03;    % k2
Ei(3) = 6.55498e+03;    % k3
Ei(4) = 1.94000e+04;    % p1
Ei(5) = 7.03000e+03;    % p2
Ei(6) = 8.10895e+04;    % r

% Frequency factor
Ai = zeros(1,6);
Ai(1) = 42.9505;        % k1
Ai(2) = 2.47851;        % k2
Ai(3) = 2.48912;        % k3
Ai(4) = 2.45614e+02;    % p1
Ai(5) = 2.61844;        % p2
Ai(6) = 1.13443e+13;    % r

Temperature = 20;       % Temperature
Ri = 8.314;             % Gas constant

%%% Parameter values for Goodwin model %%%
k1 = Ai(1)*exp(-Ei(1)/Ri/(273+Temperature));
k2 = Ai(2)*exp(-Ei(2)/Ri/(273+Temperature));
k3 = Ai(3)*exp(-Ei(3)/Ri/(273+Temperature));
p1 = Ai(4)*exp(-Ei(4)/Ri/(273+Temperature));
p2 = Ai(5)*exp(-Ei(5)/Ri/(273+Temperature));
v = Ai(6)*exp(-Ei(6)/Ri/(273+Temperature));
n = 20;
K = 0.0184;

%%% Paramter values for external force %%%
L = 0.001;      % Amplitude of external force
Omega = 0.26;   % Angular velocity of external force

%%% Initial values %%%
x0 = [0.105066740477856 0.021160486166078 0.016385059112311 0];

%%% Simulation for Goodwin model by Runge-Kutta method of order 4 %%%
[t, x] = rungekutta4(@Goodwin_Hill_externalf,tspan,x0,dt);

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

%%% Timeseries of x3 %%%
figure
plot(t,x(:,3))
xlim([0 100])
box on

fprintf('period = %4.3f\n', period);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Goodwin model with external force %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Goodwin_Hill_externaldef = Goodwin_Hill_externalf(~,x)

global k1 k2 k3 p1 p2 v n K L Omega

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

x1_dot = v/(1+(x3/K)^n)-k1*x1+L*cos(x4);
x2_dot = p1*x1-k2*x2;
x3_dot = p2*x2-k3*x3;
x4_dot = Omega;

output(1) = x1_dot;
output(2) = x2_dot;
output(3) = x3_dot;
output(4) = x4_dot;

% return a column vector
output = output(:);

Goodwin_Hill_externaldef = output;
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