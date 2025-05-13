%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Generalize Harmonic Aanalysis (GHA) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameters for GHA %%%
Fs = 1/dt;                  % Sampling frequency
Pd = period;                % Period
FramePoint = length(y);     % Signal length
ExtractNum = Enum;          % Extract number of frequency

frequency = zeros(1, ExtractNum);
amplitude = zeros(1, ExtractNum);
phs = zeros(1, ExtractNum);
sphase = zeros(1, ExtractNum);
cphase = zeros(1, ExtractNum);

xx = double(y(1:FramePoint)) ;
xx = xx';
nn = 0:1/Fs:(FramePoint-1)/Fs;

af = zeros(100, 1) ;
bf = zeros(100, 1) ;

%%% loop for extraction %%%
for j=1:ExtractNum
    %%% loop for every frequency to minimize error function %%%
    for fR=1:100
        cc = sum(xx.*cos(2*pi*fR/Pd*nn))/Fs;
        ss = sum(xx.*sin(2*pi*fR/Pd*nn))/Fs;

        %%% Amplitude %%%
        zz = 4*pi*fR/Pd*FramePoint / Fs ;
        P = (zz + sin(zz)) / (8*pi*fR/Pd) ;
        Q = (zz - sin(zz)) / (8*pi*fR/Pd) ;
        R = (1. - cos(zz)) / (8*pi*fR/Pd) ;

        aaf = (Q*cc - R*ss) / (P*Q - R*R) ;
        bbf = (P*ss - R*cc) / (P*Q - R*R) ;


        %%% Error function %%%
        res = xx - (aaf*cos(2*pi*fR/Pd*nn) + bbf*sin(2*pi*fR/Pd*nn));
        pw = sum(res.^2);

        if fR==1
            Emin = pw ;
            extf = fR ;
            exta = aaf ;
            extb = bbf ;
        else
            if Emin > pw
                Emin = pw ;
                extf = fR ;
                exta = aaf ;
                extb = bbf ;
            end
        end


    end

    %%% Residue %%%
    xx = xx - (exta*cos(2*pi*extf/Pd*nn) + extb*sin(2*pi*extf/Pd*nn));

    af(extf) = af(extf) + exta ;
    bf(extf) = bf(extf) + extb ;
end

mR = 0 ;
fk = zeros(ExtractNum,1);
for fR=1:100
    pw = af(fR)*af(fR)+bf(fR)*bf(fR) ;
    if pw>0
        mR = mR+1 ;
        fk(mR) = fR ;
    end
    if mR==ExtractNum
        break ;
    end
end

for j=1:mR
    frequency(1, j) = fk(j)/Pd ;
    amplitude(1, j) = sqrt( af(fk(j))*af(fk(j)) + bf(fk(j))*bf(fk(j)) ) ;
    phs(1, j) = atan2( bf(fk(j)), af(fk(j)) ) * 360 / (2*pi) ;
    sphase(1, j) = atan2( af(fk(j)), bf(fk(j)) );
    cphase(1, j) = atan2( bf(fk(j)), af(fk(j)) );
end

%%% Absolute value of Fourier coeffient %%%
Asquare = amplitude.^2/4;
au = Asquare.*fk'.^4;
ad = Asquare.*fk'.^2;
at = Asquare;