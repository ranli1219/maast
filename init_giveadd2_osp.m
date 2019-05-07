function init_giveadd2_osp()
%*************************************************************************
%*     Copyright c 2009 The board of trustees of the Leland Stanford     *
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

global GIVE2_CHI2_PFA GIVE2__CHI2_PMD
global GIVE2_SIG2_DECORR GIVE2_CHI2_THRESH GIVE2_CHI2_LOWER GIVE2_CHI2_NOM
global GIVE2_RIRREG2_FLOOR GIVE2_CHI2_LWR_UPM
global GIVE2_N_MIN GIVE2_N_PTS GIVE2_RMAX GIVE2_RMIN
global GIVE2_RCMX_N_MIN GIVE2_RCMX_N_PTS GIVE2_RCMX_RMAX GIVE2_RCMX_RMIN
global GIVE2_UNDERSAMP_PTS GIVE2_UPMUNSAMP_PTS
global GIVE2_KHMI_KA GIVE2_GIVEI_FLOOR
global GIVE2_KA_RATIO
global GIVE2_SIG2_L1L2REC GIVE2_SIG2_L1L2SAT
global GIVE2_MEX_SLOPE1 GIVE2_MEX_SLOPE2 GIVE2_MEX_LLAT1 GIVE2_MEX_LLAT2
global GIVE2_ANT_A1 GIVE2_ANT_B1 GIVE2_ANT_C1 GIVE2_ANT_ALPHA1 
global GIVE2_ANT_A2 GIVE2_ANT_B2 GIVE2_ANT_C2 GIVE2_ANT_ALPHA2
global RTR_FLAG RTR_MULT TRUTH_FLAG IPP_SPREAD_FLAG BRAZPARAMS


GIVE2_CHI2_PFA = 10^(-3);
GIVE2_CHI2_PMD = 10^(-6);
GIVE_CHI2_PHMI = 10^(-10);
if BRAZPARAMS
    GIVE2_SIG2_DECORR=2^2;
else
    GIVE2_SIG2_DECORR=(0.35)^2;
end;
fprintf('Sig2 Decorr = %f\n', GIVE2_SIG2_DECORR);

GIVE2_KHMI_KA=5.592;

GIVE2_CHI2_THRESH  = chi2inv(1-GIVE2_CHI2_PFA,1:200)';
%GIVE2_CHI2_LOWER   = chi2inv(GIVE2_CHI2_PMD,1:200)';
GIVE2_CHI2_LOWER   = RTR_table(200, GIVE_CHI2_PHMI, GIVE2_KHMI_KA);
GIVE2_CHI2_LWR_UPM = chi2inv(0.5,1:200)';
GIVE2_CHI2_NOM     = chi2inv(.25,1:200)';

GIVE2_RIRREG2_FLOOR = 1;

GIVE2_N_MIN=10;
GIVE2_N_PTS=30;
GIVE2_RMAX=2100000;
GIVE2_RMIN=800000;

                   
GIVE2_KA_RATIO  = GIVE2_KHMI_KA/5.33;

GIVE2_SIG2_L1L2REC = (0.229)^2;
GIVE2_SIG2_L1L2SAT = (0.192)^2;

GIVE2_GIVEI_FLOOR = 10;  %GIVE floor set to 3 m

GIVE2_MEX_LINE1 = [[21.75 13]; [-130 -80]];
GIVE2_MEX_LINE2 = [[13 13]; [-80 -60]];
GIVE2_MEX_SLOPE1 = (GIVE2_MEX_LINE1(1,1) - GIVE2_MEX_LINE1(1,2))/ ...
                   (GIVE2_MEX_LINE1(2,1) - GIVE2_MEX_LINE1(2,2));
GIVE2_MEX_SLOPE2  = (GIVE2_MEX_LINE2(1,1) - GIVE2_MEX_LINE2(1,2))/ ...
                   (GIVE2_MEX_LINE2(2,1) - GIVE2_MEX_LINE2(2,2));
GIVE2_MEX_LLAT1 = GIVE2_MEX_LINE1(1,1) - GIVE2_MEX_SLOPE1*GIVE2_MEX_LINE1(2,1);
GIVE2_MEX_LLAT2 = GIVE2_MEX_LINE2(1,1) - GIVE2_MEX_SLOPE2*GIVE2_MEX_LINE1(2,1);

GIVE2_ANT_A1     = 0.030;
GIVE2_ANT_B1     = 0.20;
GIVE2_ANT_C1     = 1.5;
GIVE2_ANT_ALPHA1 = 21*pi/180;

GIVE2_ANT_A2     = 0.030;
GIVE2_ANT_B2     = 0.2;
GIVE2_ANT_C2     = 1.5;
GIVE2_ANT_ALPHA2 = 21*pi/180;

GIVE2_UPMUNSAMP_PTS = [  [0.700 2100 3.32580];...
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
                         [0.075  800 0.16686];...
                         [0.000    0 0.0000]];
                     
                     
                     
GIVE2_UNDERSAMP_PTS = [  [0.600 2100 3.5538];...
                         [0.575 2100 3.4853];...
                         [0.550 2100 3.3971];...
                         [0.625 1750 3.3354];...
                         [0.525 2100 3.2293];...
                         [0.600 1750 3.2092];...
                         [0.575 2050 3.1136];...
                         [0.550 2050 3.0338];...
                         [0.500 2100 3.0198];...
                         [0.575 2000 2.8105];...
                         [0.550 2000 2.6592];...
                         [0.475 2100 2.5433];...
                         [0.575 1450 2.4526];...
                         [0.450 2100 1.9393];...
                         [0.500 2050 1.8283];...
                         [0.525 2000 1.7748];...
                         [0.450 2050 1.7299];...
                         [0.525 1950 1.7174];...
                         [0.475 2000 1.6768];...
                         [0.450 2000 1.6370];...
                         [0.525 1900 1.6213];...
                         [0.500 1900 1.5942];...
                         [0.425 2000 1.5385];...
                         [0.475 1950 1.5368];...
                         [0.375 2100 1.5086];...
                         [0.450 1950 1.4791];...
                         [0.425 1950 1.4546];...
                         [0.500 1800 1.4396];...
                         [0.425 1900 1.4132];...
                         [0.400 1900 1.4111];...
                         [0.400 1850 1.3649];...
                         [0.375 2050 1.3527];...
                         [0.400 1800 1.2808];...
                         [0.525 1750 1.2500];...
                         [0.575 1350 1.2283];...
                         [0.525 1350 1.1863];...
                         [0.375 1900 1.1765];...
                         [0.375 1800 1.1575];...
                         [0.500 1750 1.0942];...
                         [0.500 1500 1.0942];...
                         [0.600 1300 1.0660];...
                         [0.575 1250 1.0624];...
                         [0.550 1250 1.0430];...
                         [0.400 1700 1.0389];...
                         [0.500 1450 1.0209];...
                         [0.500 1300 1.0177];...
                         [0.350 2050 1.0086];...
                         [0.525 1250 1.0042];...
                         [0.400 1150 0.9884];...
                         [0.400 1100 0.9720];...
                         [0.375 1100 0.9157];...
                         [0.350 2000 0.9120];...
                         [0.275  950 0.7831];...
                         [0.250  950 0.7591];...
                         [0.050 2100 0.7069];...
                         [0.200  950 0.6953];...
                         [0.175  850 0.6126];...
                         [0.075 1850 0.5840];...
                         [0.075 1750 0.5380];...
                         [0.050 1850 0.4778];...
                         [0.075 1600 0.4223];...
                         [0.025 2100 0.4191];...
                         [0.125  850 0.4121];...
                         [0.050 1800 0.3732];...
                         [0.050 1750 0.3688];...
                         [0.075  900 0.3029];...
                         [0.125  800 0.2239];...
                         [0.100  800 0.1808];...
                         [0.075  850 0.1457];...
                         [0.050 1700 0.1088];...
                         [0.075  800 0.0921];...
                         [0.000    0 0.0000]];
                     