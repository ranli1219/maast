function init_udre2_osp()
%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
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
%function init_udre_osp()
% 
% assigns values used by the UDRE ADD
%
% References: 	UDRE ADD Version 8

% Created 2001 Apr03 by Wyant Chan
% Modified 2001 Apr18 by Todd Walter
% Modified 2001 Jun19 by Todd Walter
% Modified 2001 Nov01 by Todd walter
% Added Antenna bias parameters 2004 May21 by Todd Walter
% Updated to Release 3A1 2013 Aug 30 by Todd Walter

global UDRE2_SIG2_L1L2 UDRE2_SIG2_WRECLK
global UDRE2_SIG_SV_H UDRE2_SIG_SV_C UDRE2_SIG_SV_L UDRE2_SIG_SV_T ...
        UDRE2_SIG_DOT_H UDRE2_SIG_DOT_C UDRE2_SIG_DOT_L UDRE2_T_PREDEG
global UDRE2_PHMI UDRE2_KHMI UDRE2_KFA
global UDRE2_MIN_NWRS UDRE2_UDREI_FLOOR
global UDRE2_DW_ELEV UDRE2_DW1 UDRE2_DW2 UDRE2_DW3 UDRE2_MT28_MIN_CNMP2
global UDRE2_PHI_LOWER_BND UDRE2_PHI_CUTOFF
global UDRE2_CLA_GPS UDRE2_ALPHA_TROP
global UDRE2_K1 UDRE2_K2 UDRE2_K3 UDRE2_K4 ...
       UDRE2_K5 UDRE2_K6 UDRE2_K7 UDRE2_K8 
		
		
UDRE2_SIG2_L1L2 = (.192/3.29)^2;
UDRE2_SIG2_WRECLK = (0.35/3.29)^2;

UDRE2_SIG_SV_H  = 2.61;
UDRE2_SIG_SV_C  = 5.45;
UDRE2_SIG_SV_L  = 13.25;
UDRE2_SIG_SV_T  = 100;
UDRE2_SIG_DOT_H = 1.3e-3;
UDRE2_SIG_DOT_C = 0.9e-3;
UDRE2_SIG_DOT_L = 0.5e-3;
UDRE2_T_PREDEG  = 360;

UDRE2_PHMI      = 4.5e-10;    % prob of HMI allocation to UDRE2 monitor
UDRE2_KHMI      = 5.33;       % matches MOPS KV_PA
UDRE2_KFA       = 3.29/2;     % set empirically to trade avail. and continuity

UDRE2_MIN_NWRS          = 3;  % Minimum number of visible WRS
UDRE2_UDREI_FLOOR       = 5 +1;  % Minimum Allowable UDREI (matlab runs 1:16 not 0:15 as in the MOPS)

UDRE2_DW_ELEV =  7*pi/180;    % radians below which to start deweighting
UDRE2_DW1     =  0.238060528495;
UDRE2_DW2     = -1.17983911751;
UDRE2_DW3     = -1.055141155463;

UDRE2_ALPHA_TROP = (1/0.5)^2 - 1; %tropospheric adjustment for threshold

UDRE2_MT28_MIN_CNMP2 = (0.45)^2;  % minimum CNMP value for MT28 in m


UDRE2_PHI_LOWER_BND = 0.1;        % lowest input value to use in phi func.
UDRE2_PHI_CUTOFF    = 5.0;        % largest input value to use in phi func.

UDRE2_CLA_GPS = 1.1;  % critical limit adjustment
UDRE2_K1 = 1.5;       % antenna bias parameters
UDRE2_K2 = 0.03;
UDRE2_K3 = 0.20;
UDRE2_K4 = 21;
UDRE2_K5 = 1.5;
UDRE2_K6 = 0.03;
UDRE2_K7 = 0.2;
UDRE2_K8 = 21;










