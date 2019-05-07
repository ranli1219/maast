function init_give_osp()
%*************************************************************************
%*     Copyright c 2001 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     The algorithms in this script file may contain elements           *
%*     considered proprietary by Raytheon. This file shall not be        *
%*     disclosed to the public, to any vendor not identified as a        *
%*     consultant subcontractor to the FAA or the prime contractor       *
%*     (Raytheon), or to any other government (or government agency) not *
%*     participating in WAAS without the prior written permission of the *
%*     FAA WAAS Contracting Officer                                      *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
%function init_give_osp()

%modified 05 Nov 2001 by Todd Walter
%modified 24 Apr 2003 by Juan Blanch

global KRIG_CHI2_PFA KRIG__CHI2_PMD
global KRIG_SIG2_DECORR KRIG_R2_IRREG 
global KRIG_N_MIN KRIG_N_PTS KRIG_RMAX KRIG_RMIN
global KRIG_UNDERSAMP_PTS
global KRIG_SIG2_ROT0
global KRIG_KHMI_KA KRIG_GIVEI_FLOOR
global KRIG_KIGD_FAC KRIG_CEILING_P KRIG_CEILING_T KRIG_CEILING_LAT
global KRIG_SIG2_L1L2REC KRIG_SIG2_L1L2SAT
global KRIG_NUGGET_EFFECT KRIG_DIST2DOUBLE KRIG_VARIANCE
global KRIG_CHI2_THRESH;
global RTR_FLAG RTR_MULT TRUTH_FLAG IPP_SPREAD_FLAG BRAZPARAMS KRIG_ORDER

KRIG_CHI2_PFA=10^(-3);
KRIG_CHI2_PMD=(2.25*10^(-4));

KRIG_SIG2_DECORR=(0.35)^2;

KRIG_CHI2_THRESH = chi2inv(1-KRIG_CHI2_PFA,1:100)';

%Variogram parameters:
% KRIG_NUGGET_EFFECT is the intercept of the variogram at distance 0
% (can be thought as a sigma_decorr at distance 0)
% KRIG_DIST2DOUBLE is the distance over which the variogram doubles
% the nugget effect.
% KRIG_VARIANCE is the overall variance of the random field.
if BRAZPARAMS
    KRIG_NUGGET_EFFECT = (1).^2;
    KRIG_DIST2DOUBLE = 600000;
    KRIG_VARIANCE = 16;
else
    KRIG_NUGGET_EFFECT=.04;
    KRIG_DIST2DOUBLE=650000;
    KRIG_VARIANCE = 2;
end;
fprintf('Nugget = %f  Dist2Double = %f  Variance = %f\n', KRIG_NUGGET_EFFECT, KRIG_DIST2DOUBLE, KRIG_VARIANCE);

% Kriging order
KRIG_ORDER = 1;

if RTR_FLAG
    load('RTR_MULT.mat');
    fprintf('Real Time R-irreg being used\n');
else
    fprintf('Normal R-irreg being used\n');
end;

if IPP_SPREAD_FLAG
    fprintf('IPP Spread Metric Being Used\n');
else
    fprintf('Relative Centroid Metric Being Used\n');
end;

KRIG_R2_IRREG=(chi2inv(1-KRIG_CHI2_PFA,1:200)./chi2inv(KRIG_CHI2_PMD,1:200))';

KRIG_N_MIN=10;
KRIG_N_PTS=30;
KRIG_RMAX=2100000;
KRIG_RMIN=800000;

if (IPP_SPREAD_FLAG)
% Latest version of spread metric as of 8/15/03
% A = 7, R not used, Overlap = 2, get rid of points in threat domain
% radius = 800, nug = .04 dist2dble = .65
    KRIG_UNDERSAMP_PTS = [[0.675 2100 5.8889];...
                          [0.625 2100 5.3812];...
                          [0.600 2100 3.5949];...
                          [0.550 2100 3.5694];...
                          [0.500 2100 2.5087];...
                          [0.500 1950 2.4166];...
                          [0.500 1800 2.3571];...
                          [0.475 2100 2.1252];...
                          [0.500 1700 2.1030];...
                          [0.475 1650 2.0593];...
                          [0.425 1550 1.5400];...
                          [0.375 1450 1.4948];...
                          [0.425 1200 0.4521];...
                          [0.400 1150 0.2416];...
                          [0.000    0 0.0000]];
else
 %1250 km radius nug=.04, dist2dble=650km    
 KRIG_UNDERSAMP_PTS=[0            0            0
       0.375         2100       0.3044
         0.4         2100       0.4931
       0.425         2100       1.4854
        0.45         2100       2.1125
        0.45         1750       1.2673
       0.475         2100       2.2655
       0.525         1450       1.5151
       0.525         1250       0.3658
        0.55         2100       3.6102
        0.55         1250       0.4731
       0.575         1700       2.1266
       0.575         1650       2.0825
       0.575         1550       1.5617
         0.6         1950       2.4429
         0.6         1800       2.3829
        0.65         2100        5.442
         0.8         2100       5.9835];
 end;

% Test 800 km radius, use ipp spread metric (full tm run)
% radius = 800 or 1200
%       KRIG_UNDERSAMP_PTS = [[  0.700 2100 5.8889];
%                               [0.675 2100 5.3812];
%                               [0.650 2100 3.5949];
%                               [0.600 2100 3.5694];
%                               [0.550 2100 2.5087];
%                               [0.550 1950 2.4166];
%                               [0.525 2100 2.3597];
%                               [0.525 1800 2.3571];
%                               [0.525 1700 2.1030];
%                               [0.525 1650 2.0593];
%                               [0.500 2100 1.8715];
%                               [0.500 1650 1.5419];
%                               [0.450 1550 1.5400];
%                               [0.425 1450 1.4948];
%                               [0.475 1200 0.4521];
%                               [0.450 1150 0.2416];
%                               [0.000    0 0.0000]];
%                         
%1200 km radius nug=.05, dist2dble=1000km
%KRIG_UNDERSAMP_PTS=[[0  800    0];
%          [0.3         2100       0.7973/1.1];
%        [0.375         2000       0.5437/1.1];
%        [  0.4         2100       0.9609/1.1];
%        [0.425         2100       1.9855/1.1];
%        [ 0.45         2100       2.3995/1.1];
%        [ 0.45         1800       0.4002/1.1];
%        [ 0.45         1650       0.1727/1.1];
%       [ 0.475         2100       2.9299/1.1];
%       [ 0.475         1800       0.5315/1.1];
%       [ 0.475         1500       0.4711/1.1];
%       [   0.5         1700       0.6545/1.1];
%       [   0.5         1400       0.5866/1.1];
%       [ 0.525         1500       1.6395/1.1];
%       [ 0.525         1450       1.6314/1.1];
%      [ 0.525         1050       0.5869/1.1];
%      [  0.55         2100       3.7526/1.1];
%       [ 0.575         1550       1.6994/1.1];
%       [ 0.575         1200       0.7117/1.1];
%       [ 0.675         2100       4.9066/1.1];
%       [   0.7         2000       1.8012/1.1];
%       [ 0.725         2100       5.8523/1.1];
%       [ 0.725         1900       1.7221/1.1];
%       [ 0.725         1850       1.7178/1.1];
%       [  0.75         1850       1.7389/1.1]];

               
   
%1200 km radius nug=.05, dist2dble=1000km 
 %KRIG_UNDERSAMP_PTS=[          0          800            0;
 %        0.35         2100        0.649;
 %       0.425         2100       1.4926;
 %        0.45         2100       2.1086;
 %        0.45         1750       1.2552;
 %       0.475         2100        2.268;
 %         0.5         1400       1.4965;
 %       0.525         1500         1.51;
 %       0.525         1450        1.509;
 %        0.55         2100       3.6145;
 %        0.55         1850       2.2333;
 %        0.55         1150        0.277;
 %       0.575         2100       3.6394;
 %       0.575         1700       2.1159;
 %       0.575         1650       2.0724;
 %       0.575         1550       1.5604;
 %       0.575         1200       0.4297;
 %         0.6         1800        2.389;
 %       0.625         2100       5.1967;
 %       0.625         1300       0.4686;
 %        0.65         2100       5.4579;
 %         0.8         2100       6.0257];
      

%kriging metric 800 km .05 dist2double=1000km
%KRIG_UNDERSAMP_PTS=[    0            800           0;
%        0.259         2100       1.2823;
%        0.259         1750       1.2552;
%        0.273         1500       0.4746;
%        0.273         1150       0.2063;
%        0.287         1450        1.509;
%        0.287         1400       1.4965;
%        0.294         2000       2.0796;
%        0.294         1800       2.0649;
%        0.294         1750       2.0309;
%        0.294         1700       1.9132;
%        0.294         1500         1.51;
%        0.294         1200       0.4297;
%        0.301         1950       2.5645;
%        0.301         1800        2.389;
%        0.301         1700       2.1159;
%        0.301         1650       2.0724;
%        0.301         1550       1.5604;
%        0.301         1300       0.4686;
%        0.399         2100         3.53;
%        0.406         2100       3.6145;
%        0.476         2100       5.4579;
%        0.686         2100       6.0257];
%KRIG_UNDERSAMP_PTS=[    0          800            0
%        0.273         2100       0.3042
%         0.28         2100       0.5801
%         0.28         1850       0.3213
%        0.287         1750       1.0383
%        0.294         1750       1.2552
%        0.301         1850       1.5955
%        0.308         1750       1.8321
%        0.308         1600       0.7135
%        0.308         1500       0.4746
%        0.308         1150       0.2063
%        0.315         1750       2.0309
%        0.315         1700       1.9132
%        0.315         1550       1.0821
%        0.322         1850       2.2333
%        0.322         1650       2.0724
%        0.322         1550         1.16
%        0.322         1200       0.4297
%        0.329         1950       2.4503
%        0.329         1800        2.389
%        0.329         1700       2.1159
%        0.329         1300       0.4686
%        0.336         1950       2.5645
%        0.336         1500       0.5084
%        0.343         1500       1.0473
%         0.35         1450        1.509
%         0.35         1400       1.4965
%        0.357         1500         1.51
%        0.364         1600        1.548
%        0.378         1550       1.5604
%        0.399         2100       2.7034
%        0.406         2100       2.7643
%        0.413         2100       3.3405
%         0.42         2100         3.53
%        0.427         2100       3.6145
%        0.532         2100       5.4579
%        0.686         2100       6.0257];

%1500 km radius nug=.04, dist2dble=650km
 %KRIG_UNDERSAMP_PTS=[          0            0            0
 %                          0.375         2100       0.3044
 %         0.4         2100       0.4931
 %       0.425         2100       1.4854
 %        0.45         2100       2.1125
 %        0.45         1750       1.2673
 %       0.475         2100       2.2655
 %         0.5         1500       1.5151
 %        0.55         2100       3.6102
 %       0.575         1700       2.1266
 %       0.575         1650       2.0825
 %       0.575         1550       1.5617
 %         0.6         1950       2.4429
 %         0.6         1800       2.3829
 %        0.65         2100        5.442
 %         0.8         2100       5.9835];



 %kriging metric- 1300 km minimum radius nug=.04, dist2double=650
 %    KRIG_UNDERSAMP_PTS=[       0          800            0
 %       0.252         2100       0.4456
 %       0.259         2100       0.4854
 %       0.266         2100       0.6612
 %       0.266         1850       0.3978
 %       0.273         1750       1.0558
 %        0.28         1750       1.2673
 %       0.287         2100        1.291
 %       0.294         1850       1.5965
 %       0.294         1650       0.6851
 %       0.294         1550        0.371
 %       0.301         1750       1.8357
 %       0.301         1600        0.739
 %       0.301         1550       0.5581
 %       0.301         1300        0.399
 %       0.308         1750       1.9633
 %       0.308         1700       1.9216
 %       0.308         1550        1.089
 %       0.315         1850        2.197
 %       0.315         1650       2.0825
 %       0.315         1300       0.5095
 %       0.322         1950       2.4429
 %       0.322         1800       2.3829
 %       0.322         1700       2.1266
 %       0.329         2100       2.4873
 %       0.336         1500       0.6376
 %       0.343         1450       1.5151
 %       0.364         1600         1.55
 %       0.378         1550       1.5617
 %       0.399         2100       2.7013
 %       0.413         2100       3.3369
 %        0.42         2100       3.5237
 %       0.434         2100       3.6102
 %       0.546         2100        5.442
 %       0.686         2100       5.9835];


KRIG_SIG2_ROT0=(1.62/5.33)^2;

KRIG_KHMI_KA=5.592/5.33;

KRIG_KIGD_FAC=1/5.33;
KRIG_CEILING_P=[[24.3 9.72]; [24.3 24.3]; [72.9 9.72]; [24.3 9.72]]';
KRIG_CEILING_T=[[6 7]; [10 22]; [23 23]; [24 24]]'*3600;
KRIG_CEILING_LAT=50;

KRIG_SIG2_L1L2REC = (0.229)^2;
KRIG_SIG2_L1L2SAT = (0.192)^2;

KRIG_GIVEI_FLOOR = 10;  %GIVE floor set to 3 m
%KRIG_GIVEI_FLOOR = 0;  %testing with GIVE floor set to 2.1 m


