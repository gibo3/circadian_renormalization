%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lotka-Volterra eauation simulation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

datetime('now')

global a b eps epsy

%%% Parameters for simulation %%%
tspan = [0 8];      % timespan
dt = 0.0001;        % time step size                                
Enum = 50;          % extract number of frequency

%%% Parameter values %%%
a = 10.3581;
b = 5.1160;
eps = 11.7723;
epsy = 1.0192;

%%% Initial values %%%
x0 = [0.606900218072858 0.879881096491837];

%%% Simulation for Lotka-Volterra model by Runge-Kutta method of order 4 %%%
[t, x] = rungekutta4(@LotkaVolterramodelf,tspan,x0,dt);

%%% Find peak index %%%
Apeakx = zeros(40000,1);
ipeakx = zeros(40000,1);

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
NS = sqrt(sum(ad(1:5))/sum(at(1:5)));

%%% Timeseries of x %%%
figure
plot(t,x(:,1))
xlim([0 3])
box on

fprintf('period = %4.3f, NS = %4.4f\n', period, NS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lotka-Volterra Model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LotkaVolterramodeldef = LotkaVolterramodelf(~,x)

global a b eps epsy

X = x(1);
Y = x(2);

X_dot = a*X-eps*X*Y;
Y_dot = -b*Y+epsy*X*Y;

output(1) = X_dot;
output(2) = Y_dot;

% return a column vector
output = output(:);

LotkaVolterramodeldef = output;
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