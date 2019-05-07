function init_giveadd1_osp()
%*************************************************************************
%*     Copyright c 2007 The board of trustees of the Leland Stanford     *
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
%function init_giveadd1_osp()

%created 18 May 2004 by Todd Walter
%updated 29 June 2006 by Todd Walter to align it to the ADD

global GIVE1_CHI2_PFA GIVE1__CHI2_PMD
global GIVE1_SIG2_DECORR GIVE1_CHI2_THRESH GIVE1_CHI2_LOWER GIVE1_CHI2_NOM
global GIVE1_RIRREG2_FLOOR GIVE1_CHI2_LWR_UPM
global GIVE1_N_MIN GIVE1_N_PTS GIVE1_RMAX GIVE1_RMIN
global GIVE1_UNDERSAMP_PTS GIVE1_UPMUNSAMP_PTS
global GIVE1_SIG2_ROT0
global GIVE1_KHMI_KA GIVE1_GIVEI_FLOOR
global GIVE1_KA_RATIO
global GIVE1_SIG2_L1L2REC GIVE1_SIG2_L1L2SAT
global GIVE1_MEX_SLOPE1 GIVE1_MEX_SLOPE2 GIVE1_MEX_LLAT1 GIVE1_MEX_LLAT2
global GIVE1_ANT_A1 GIVE1_ANT_B1 GIVE1_ANT_C1 GIVE1_ANT_ALPHA1 
global GIVE1_ANT_A2 GIVE1_ANT_B2 GIVE1_ANT_C2 GIVE1_ANT_ALPHA2
global RTR_FLAG RTR_MULT TRUTH_FLAG IPP_SPREAD_FLAG BRAZPARAMS


GIVE1_CHI2_PFA = 10^(-3);
GIVE1_CHI2_PMD = 10^(-6);
GIVE_CHI2_PHMI = 10^(-10);
if BRAZPARAMS
    GIVE1_SIG2_DECORR=2^2;
else
    GIVE1_SIG2_DECORR=(0.35)^2;
end;
fprintf('Sig2 Decorr = %f\n', GIVE1_SIG2_DECORR);

GIVE1_KHMI_KA=5.592;

GIVE1_CHI2_THRESH  = chi2inv(1-GIVE1_CHI2_PFA,1:200)';
%GIVE1_CHI2_LOWER   = chi2inv(GIVE1_CHI2_PMD,1:200)';
GIVE1_CHI2_LOWER   = RTR_table(200, GIVE_CHI2_PHMI, GIVE1_KHMI_KA);
GIVE1_CHI2_LWR_UPM = chi2inv(0.5,1:200)';
GIVE1_CHI2_NOM     = chi2inv(.25,1:200)';

GIVE1_RIRREG2_FLOOR = 1;

GIVE1_N_MIN=10;
GIVE1_N_PTS=30;
GIVE1_RMAX=2100000;
GIVE1_RMIN=800000;


fprintf('Relative Centroid Metric being used\n');
% rcm, shell height = 350
% From Threat model with inclusion of ADD1 and October & November storms
                     
GIVE1_UNDERSAMP_PTS = [  [0.725 2100 3.4463];...
                         [0.675 2100 3.3882];...
                         [0.675 2000 2.9979];...
                         [0.650 2000 2.6964];...
                         [0.675 1950 2.5717];...
                         [0.675 1900 2.5316];...
                         [0.575 1450 2.4526];...
                         [0.575 1350 1.2283];...
                         [0.525 1350 1.1863];...
                         [0.600 1300 1.0660];...
                         [0.575 1250 1.0624];...
                         [0.550 1250 1.0430];...
                         [0.500 1300 1.0177];...
                         [0.525 1250 1.0042];...
                         [0.400 1150 0.9889];...
                         [0.400 1100 0.9768];...
                         [0.375 1100 0.9468];...
                         [0.275 1900 0.8927];...
                         [0.275 1850 0.8750];...
                         [0.250  950 0.8452];...
                         [0.225  950 0.7476];...
                         [0.075 2100 0.7309];...
                         [0.050 2100 0.7069];...
                         [0.150  850 0.6626];...
                         [0.100 1850 0.5436];...
                         [0.075 1850 0.5206];...
                         [0.125  850 0.5097];...
                         [0.075 1750 0.4985];...
                         [0.150  800 0.4620];...
                         [0.050 1750 0.4545];...
                         [0.050 1700 0.4401];...
                         [0.025 2100 0.4191];...
                         [0.075 1600 0.3708];...
                         [0.100  950 0.2679];...
                         [0.075 1350 0.2204];...
                         [0.075  850 0.1769];...
                         [0.000    0 0.0000]];
                     
GIVE1_SIG2_ROT0 = 0.01;

GIVE1_KA_RATIO  = GIVE1_KHMI_KA/5.33;

GIVE1_SIG2_L1L2REC = (0.229)^2;
GIVE1_SIG2_L1L2SAT = (0.192)^2;

GIVE1_GIVEI_FLOOR = 10;  %GIVE floor set to 3 m

GIVE1_MEX_LINE1 = [[26.75 18]; [-130 -80]];
GIVE1_MEX_LINE2 = [[18 18]; [-80 -60]];
GIVE1_MEX_SLOPE1 = (GIVE1_MEX_LINE1(1,1) - GIVE1_MEX_LINE1(1,2))/ ...
                   (GIVE1_MEX_LINE1(2,1) - GIVE1_MEX_LINE1(2,2));
GIVE1_MEX_SLOPE2  = (GIVE1_MEX_LINE2(1,1) - GIVE1_MEX_LINE2(1,2))/ ...
                   (GIVE1_MEX_LINE2(2,1) - GIVE1_MEX_LINE2(2,2));
GIVE1_MEX_LLAT1 = GIVE1_MEX_LINE1(1,1) - GIVE1_MEX_SLOPE1*GIVE1_MEX_LINE1(2,1);
GIVE1_MEX_LLAT2 = GIVE1_MEX_LINE2(1,1) - GIVE1_MEX_SLOPE2*GIVE1_MEX_LINE1(2,1);

GIVE1_ANT_A1     = 0.030;
GIVE1_ANT_B1     = 0.20;
GIVE1_ANT_C1     = 1.5;
GIVE1_ANT_ALPHA1 = 21*pi/180;

GIVE1_ANT_A2     = -0.015;
GIVE1_ANT_B2     = 0.36;
GIVE1_ANT_C2     = 1.3;
GIVE1_ANT_ALPHA2 = 47*pi/180;

GIVE1_UPMUNSAMP_PTS = [  [0.700 2100 3.32580];...
                         [0.675 2100 3.32020];...
                         [0.650 2100 3.11850];...
                         [0.625 2100 2.97880];...
                         [0.600 2100 1.62790];...
                         [0.575 2100 1.57250];...
                         [0.550 2100 1.50210];...
                         [0.600 1950 1.59700];...
                         [0.575 1950 1.52370];...
                         [0.650 1900 1.40550];...
                         [0.625 1900 1.23220];...
                         [0.600 1900 1.19300];...
                         [0.050 1900 0.67717];...
                         [0.650 1850 1.00710];...
                         [0.650 1800 0.82517];...
                         [0.175 1700 0.70292];...
                         [0.075 1550 0.70010];...
                         [0.025 1550 0.29100];...
                         [0.050 1500 0.67000];...
                         [0.075 1400 0.69925];...
                         [0.050 1350 0.41194];...
                         [0.525 1300 0.74956];...
                         [0.200 1050 0.65461];...
                         [0.225 1000 0.65594];...
                         [0.175  800 0.48071];...
                         [0.125  800 0.25913];...
                         [0.075  800 0.16686]];

                     