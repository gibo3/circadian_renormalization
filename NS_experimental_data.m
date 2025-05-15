%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Genaralized Harmonics Analysis of expreimental time-series %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

%%% read time-series data %%%
C30 = csvread('timeseries_example.csv');

%%% number of data sets %%%
N = 100;
periodC30 = zeros(1,N);
NSC30 = zeros(1,N);

%%% delete overlapping time %%%
difC30 = C30(2:end,1)-C30(1:end-1,1);
istC30 = find(difC30 == 0);
uC30 = C30;
uC30(istC30,:) = [];

%%% spline interpolation %%%
uC30x = (0:1:24)';
spuC30 = spline(uC30(:,1),uC30(:,2),uC30x);

%%% normalization %%%
spnC30 = spuC30/max(abs(spuC30));

%%% noise adding %%%
noisyC30 = zeros(length(spnC30),N);

for bw = 1:N
    noisyC30(:,bw) = spnC30+0.1*(rand(length(spnC30),1)-1/2);
end

%%% spline interpolation %%%
fcC30x = (0:0.1:24)';
spnoisyC30 = zeros(length(fcC30x),N);

figure
hold on
for bd = 1:N
    spnoisyC30(:,bd) = spline(uC30x,noisyC30(:,bd),fcC30x);

    plot(fcC30x,spnoisyC30(:,bd),'c--')
end
plot(uC30x,spnC30,'b-')
hold off
box on

spleC30 = length(spnoisyC30(:,1));

%%% Parameter for GHA %%%
Enum = 100;
dt = fcC30x(2)-fcC30x(1);

%%% NS calculation %%%
for bq = 1:N
    %%% peak detection %%%
    pnC30 = zeros(1,125);
    inC30 = zeros(1,125);

    ipnC30 = 1;
    for zi = 2:spleC30-1
        if spnoisyC30(zi,bq) < spnoisyC30(zi-1,bq) && spnoisyC30(zi,bq) < spnoisyC30(zi+1,bq) && spnoisyC30(zi,bq) < -0.75  % 前後の値より大きい場合はピーク
            pnC30(ipnC30) = spnoisyC30(zi,bq);
            inC30(ipnC30) = zi;
            ipnC30 = ipnC30+1;
        end
    end

    zeroC30 = find(pnC30 == 0);
    pnC30(zeroC30) = [];
    inC30(zeroC30) = [];

    f_inC30 = find(inC30 < 100);
    l_inC30 = find(inC30 > 100);

    if length(f_inC30) == 1
        inC30_1 = inC30(f_inC30);
        pnC30_1 = pnC30(f_inC30);
    elseif isempty(f_inC30)
        inC30_1 = 1;
        pnC30_1 = spnoisyC30(1);
    else
        [pnC30_1, ka] = max(pnC30(f_inC30));
        inC30_1 = inC30(ka);
    end

    if length(l_inC30) == 1
        inC30_2 = inC30(l_inC30);
        pnC30_2 = pnC30(l_inC30);
    elseif isempty(l_inC30)
        inC30_2 = spleC30;
        pnC30_2 = spnoisyC30(end,bq);
    else
        [pnC30_2, kb] = max(pnC30(l_inC30));
        inC30_2 = inC30(l_inC30(kb));
    end

    %%% amplitude detrending %%%
    lognC30 = log([-pnC30_1 -pnC30_2]);
    rnC30 = polyfit(fcC30x([inC30_1 inC30_2]),lognC30,1);
    dnC30x = fcC30x(inC30_1:inC30_2);
    dnC30 = exp(-rnC30(1)*dnC30x).*spnoisyC30(inC30_1:inC30_2,bq);

    %%% Generalized Harmonic Analysis (GHA) %%%
    y = dnC30;
    ymean = mean(y);
    y = y-ymean;
    period = dnC30x(end)-dnC30x(1);
    periodC30(bq) = period;
    run('GHA.m')
    NSC30(bq) = sqrt(sum(au(1:3))/sum(ad(1:3)));
end

mNSC30 = mean(NSC30);

figure
hold on
swarmchart(ones(1,N),NSC30,'k*','XJitterWidth',0.6)
plot([0.7 1.3],[mNSC30 mNSC30],'k-')
xlim([0.6 1.4])
set(gca, 'XTick', []);
xlabel('example')
ylabel('NS distribution')
hold off
box on
