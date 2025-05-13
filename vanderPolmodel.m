%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% van der Pol model simulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

global eps

%%% Parameters for simulation %%%
tspan = [0 50];     % timespan
dt = 0.0001;        % time step size 
Enum = 50;          % extract number of frequency

%%% Parameter values %%%
eps = 3;

%%% Initial values %%%
x0 = [-2.023304140220708 -7.691256674572525e-05];

%%% Simulation for van der Pol model by Runge-Kutta method of order 4 %%%
[t, x] = rungekutta4(@vanderPolmodelf,tspan,x0,dt);

%%% Find peak index %%%
Apeakx = zeros(25000,1);
ipeakx = zeros(25000,1);

dipx1 = 1;
for zi = 2:length(x(:,1))-1
    if x(zi,1) > x(zi-1,1) && x(zi,1) > x(zi+1,1)
        Apeakx(dipx1) = x(zi,1);
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
y = x(:,1)-mean(x(tm1:tm2-1,1));
run('GHA.m')

%%% NS value of x %%%
NS = sqrt(sum(ad)/sum(at));

%%% Timeseries of x %%%
figure
plot(t,x(:,1))
xlim([0 20])
box on

fprintf('period = %4.3f, NS = %4.4f\n', period, NS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% van der Pol equation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vanderdef = vanderPolmodelf(~,x)

global eps

X = x(1);
Y = x(2);

X_dot = Y;
Y_dot = -X-eps*(X^2-1)*Y;

output(1) = X_dot;
output(2) = Y_dot;

% return a column vector
output = output(:);

vanderdef = output;
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
