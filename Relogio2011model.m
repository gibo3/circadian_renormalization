%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relogio et al. (2011) model calculation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

global dx1 dx2 dx3 dx5 dx6 dx7 dy1 dy2 dy3 dy4 dy5 dz1 dz2 dz3 dz4 dz5 dz6 dz7 dz8 kfx1 kdx1 kfz4 kdz4 kfz5 kdz5 kphz2 kdphz3 V1max V2max V3max V4max V5max kt1 ki1 kt2 ki2 ki21 kt3 ki3 kt4 ki4 kt5 ki5 a d g h i kp1 kp2 kp3 kp4 kp5 kiz4 kiz5 kiz6 kiz7 kiz8 kex2 kex3 b c e f f1 v w p q n m y10 y20 y30 y40 y50

%%% Parameters for simulation %%%
tspan = [0 200];        % timespan
dt = 0.01;              % time step size
Enum = 20;              % extract number of frequency

%%% Parameter values %%%
% Degradation rates for nuclear proteins or nuclear prontein complexes [/hour]
dx1 = 0.08;     % CLOCK/BMAL
dx2 = 0.06;     % PERN*/CRYN
dx3 = 0.09;     % PERN/CRYN
dx5 = 0.17;     % REV-ERBN
dx6 = 0.12;     % RORN
dx7 = 0.15;     % BMALN
% Degradation rates for mRNAs [/hour]
dy1 = 0.3;      % Per
dy2 = 0.2;      % Cry
dy3 = 2;        % Rev-Erb
dy4 = 0.2;      % Ror
dy5 = 1.6;      % Bmal
% Degradation rates for cytoplasmic proteins [/hour]
dz1 = 0.23;     % CRYc
dz2 = 0.25;     % PERc
dz3 = 0.6;      % PERc*
dz4 = 0.2;      % PERc*/CRYc
dz5 = 0.2;      % PERc/CRYc
dz6 = 0.31;     % REV-ERBc
dz7 = 0.3;      % RORc
dz8 = 0.73;     % BMALc
% Reaction rates for complex formation/dissociation
kfx1 = 2.3;     % CLOCK/BMAL-complex formation [/hour]
kdx1 = 0.01;    % CLOCK/BMAL-complex dissociation [/hour]
kfz4 = 1;       % PERc*/CRYc-complex formation [/(a.u.Ehour)]
kdz4 = 1;       % PERc*/CRYc-complex dissociation [/hour]
kfz5 = 1;       % PERc/CRYc-complex formation [/(a.u.Ehour)]
kdz5 = 1;       % PERc/CRYc-complex dissociation [/hour]
% Phosphorylation/dephosphorylation reaction rates [/hour]
kphz2 = 2;      % PERc-phosphorylation rate
kdphz3 = 0.05;  % PERc*-dephosphorylation rate
% Transcription rates [a.u./hour]
V1max = 1;      % Per
V2max = 2.92;   % Cry
V3max = 1.9;    % Rev-Erb
V4max = 10.9;   % Ror
V5max = 1;      % Bmal
% Activation/inhibition rates [a.u.]
kt1 = 3;        % Per-activation rate
ki1 = 0.9;      % Per-inhibition rate
kt2 = 2.4;      % Cry-activation rate
ki2 = 0.7;      % Cry-inhibition rate
ki21 = 5.2;     % Cry-inhibition rate
kt3 = 2.07;     % Rev-Erb-activation rate
ki3 = 3.3;      % Rev-Erb-inhibition rate
kt4 = 0.9;      % Ror-activation rate
ki4 = 0.4;      % Ror-inhibition rate
kt5 = 8.35;     % Bmal-activation rate
ki5 = 1.94;     % Bmal-inhibition rate
% Transcription fold activation (dimensionless)
a = 12;         % Per
d = 12;         % Cry
g = 5;          % Rev-Erb
h = 5;          % Ror
i = 12;         % Bmal
% Production rates [/hour]
kp1 = 0.4;      % PERc
kp2 = 0.26;     % CRYc
kp3 = 0.37;     % REV-ERBc
kp4 = 0.76;     % RORc
kp5 = 1.21;     % BMALc
% Import/Export rates [/hour]
kiz4 = 0.2;     % PERc*/CRYc
kiz5 = 0.1;     % PERc/CRYc
kiz6 = 0.5;     % REV-ERBc
kiz7 = 0.1;     % RORc
kiz8 = 0.1;     % BMALc
kex2 = 0.02;    % PERN*/CRYN
kex3 = 0.02;    % PERN/CRYN
% Hill coefficients of transcription (dimensionless)
b = 5;          % Per-activation
c = 7;          % Per-inhibition
e = 6;          % Cry-activation rate
f = 4;          % Cry-inhibition
f1 = 1;         % Cry-inhibition
v = 6;          % Rev-Erb-activation
w = 2;          % Rev-Erb-inhibition
p = 6;          % Ror-activation
q = 3;          % Ror-inhibition
n = 2;          % Bmal-activation
m = 5;          % Bmal-inhibition
% Exogenous RNA [a.u.]
y10 = 0;        % Per
y20 = 0;        % Cry
y30 = 0;        % Rev-Erb
y40 = 0;        % Ror
y50 = 0;        % Bmal

%%% Initial values %%%
x0 = [1.3655 1.2186 4.2350 0.5806 8.6076 2.3512 6.3624 1.0668 1.2639 0.0534 1.6016 6.0260 4.3823 0.2089 0.2324 0.7441 0.7217 1.3179 0.4756];

%%% Simulation for Relogio model by Runge-Kutta method of order 4 %%%
[t, x] = rungekutta4(@Relogio2011modelf,tspan,x0,dt);

%%% Find peak index %%%
Apeakx = zeros(20000,1);
ipeakx = zeros(20000,1);

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

%%% NS value of CLOCK/BMAL %%%
NS = sqrt(sum(au)/sum(ad));

%%% Timeseries of CLOCK/BMAL %%%
figure
plot(t,x(:,1))
xlim([0 100])
box on

fprintf('period = %4.2f, NS = %4.4f\n', period, NS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Relogio et al. (2011) model %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Relogio2011modeldef = Relogio2011modelf(~,x)
global dx1 dx2 dx3 dx5 dx6 dx7 dy1 dy2 dy3 dy4 dy5 dz1 dz2 dz3 dz4 dz5 dz6 dz7 dz8 kfx1 kdx1 kfz4 kdz4 kfz5 kdz5 kphz2 kdphz3 V1max V2max V3max V4max V5max kt1 ki1 kt2 ki2 ki21 kt3 ki3 kt4 ki4 kt5 ki5 a d g h i kp1 kp2 kp3 kp4 kp5 kiz4 kiz5 kiz6 kiz7 kiz8 kex2 kex3 b c e f f1 v w p q n m y10 y20 y30 y40 y50

x1 = x(1);  % CLOCK/BMAL
y3 = x(2);  % Rev-Erb
y4 = x(3);  % Ror
z6 = x(4);  % REV-ERBc
z7 = x(5);  % RORc
x5 = x(6);  % REV-ERBN
x6 = x(7);  % RORN
y5 = x(8);  % Bmal
z8 = x(9);  % BMALc
x7 = x(10); % BMALN
y1 = x(11); % Per
y2 = x(12); % Cry
z1 = x(13); % CRYc
z2 = x(14); % PERc
z3 = x(15); % PERc*
z4 = x(16); % PERc*/CRYc
z5 = x(17); % PERc/CRYc
x2 = x(18); % PERN*/CRYN
x3 = x(19); % PERN/CRYN

x1_dot = kfx1*x7 - kdx1*x1 - dx1*x1;
x2_dot = V3max*(1+g*(x1/kt3)^v)/(1+((x2+x3)/ki3)^w*(x1/kt3)^v+(x1/kt3)^v) - dy3*y3;
x3_dot = V4max*(1+h*(x1/kt4)^p)/(1+((x2+x3)/ki4)^q*(x1/kt4)^p+(x1/kt4)^p) - dy4*y4;
x4_dot = kp3*(y3+y30) - kiz6*z6 - dz6*z6;
x5_dot = kp4*(y4+y40) - kiz7*z7 - dz7*z7;
x6_dot = kiz6*z6 - dx5*x5;
x7_dot = kiz7*z7 - dx6*x6;
x8_dot = V5max*(1+i*(x6/kt5)^n)/(1+(x5/ki5)^m+(x6/kt5)^n) - dy5*y5;
x9_dot = kp5*(y5+y50) - kiz8*z8 - dz8*z8;
x10_dot = kiz8*z8 + kdx1*x1 - kfx1*x7 - dx7*x7;
x11_dot = V1max*(1+a*(x1/kt1)^b)/(1+((x2+x3)/ki1)^c*(x1/kt1)^b+(x1/kt1)^b) - dy1*y1;
x12_dot = V2max*(1+d*(x1/kt2)^e)/(1+((x2+x3)/ki2)^f*(x1/kt2)^e+(x1/kt2)^e)/(1+(x5/ki21)^f1) - dy2*y2;
x13_dot = kp2*(y2+y20) + kdz4*z4 + kdz5*z5 - kfz5*z1*z2 - kfz4*z1*z3 - dz1*z1;
x14_dot = kp1*(y1+y10) + kdz5*z5 + kdphz3*z3 - kfz5*z2*z1 - kphz2*z2 - dz2*z2;
x15_dot = kphz2*z2 + kdz4*z4 - kdphz3*z3 - kfz4*z3*z1 - dz3*z3;
x16_dot = kfz4*z1*z3 + kex2*x2 - kiz4*z4 - kdz4*z4 - dz4*z4;
x17_dot = kfz5*z1*z2 + kex3*x3 - kiz5*z5 - kdz5*z5 - dz5*z5;
x18_dot = kiz4*z4 - kex2*x2 - dx2*x2;
x19_dot = kiz5*z5 - kex3*x3 - dx3*x3;

output(1) = x1_dot;
output(2) = x2_dot;
output(3) = x3_dot;
output(4) = x4_dot;
output(5) = x5_dot;
output(6) = x6_dot;
output(7) = x7_dot;
output(8) = x8_dot;
output(9) = x9_dot;
output(10) = x10_dot;
output(11) = x11_dot;
output(12) = x12_dot;
output(13) = x13_dot;
output(14) = x14_dot;
output(15) = x15_dot;
output(16) = x16_dot;
output(17) = x17_dot;
output(18) = x18_dot;
output(19) = x19_dot;

% return a column vector
output = output(:);

Relogio2011modeldef = output;
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