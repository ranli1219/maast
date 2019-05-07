function init_udre_osp()
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

global UDRE_SIG2_SV19 UDRE_SIG2_L1L2 UDRE_SIG2_OBAD UDRE_SIG2_WRECLK
global UDRE_SIG_SV_H UDRE_SIG_SV_C UDRE_SIG_SV_L UDRE_SIG_SV_T ...
        UDRE_SIG_DOT_H UDRE_SIG_DOT_C UDRE_SIG_DOT_L UDRE_T_PREDEG
global UDRE_PHMI UDRE_KHMI UDRE_KFA
global UDRE_MIN_NWRS UDRE_MIN_NWRSCONUS UDREI_FLOOR_NOMINCONUS UDREI_FLOOR
global UDRE_DW_ELEV UDRE_DW1 UDRE_DW2 UDRE_DW3 UDRE_MT28_MIN_CNMP2
global UDRE_PHI_MU UDRE_PHI_SIG UDRE_PHI_MIN_IN 
global UDRE_CLA_GPS 
global UDRE_K1 UDRE_K2 UDRE_K3 UDRE_K4 UDRE_K5 UDRE_K6 UDRE_K7 UDRE_K8 
		
		
UDRE_SIG2_SV19 = 0;
UDRE_SIG2_L1L2 = (.192)^2;
UDRE_SIG2_OBAD = 0;
UDRE_SIG2_WRECLK = (0.35)^2;

UDRE_SIG_SV_H  = 2.61;
UDRE_SIG_SV_C  = 5.45;
UDRE_SIG_SV_L  = 13.25;
UDRE_SIG_SV_T  = 100;
UDRE_SIG_DOT_H = 1.3e-3;
UDRE_SIG_DOT_C = 0.9e-3;
UDRE_SIG_DOT_L = 0.5e-3;
UDRE_T_PREDEG  = 360;

UDRE_PHMI      = 4.5e-10;    % prob of HMI allocation to UDRE monitor
UDRE_KHMI      = 5.33;       % matches MOPS KV_PA
UDRE_KFA       = 3.29/2;     % set empirically to trade avail. and continuity

UDRE_MIN_NWRS          = 3;  % Minimum number of visible WRS
UDRE_MIN_NWRSCONUS     = 1;  % Minimum number of visible CONUS WRS
UDREI_FLOOR_NOMINCONUS = 5 +1; % Udrei floor for enough wrs but not 
                                 %  enough conus wrs
UDREI_FLOOR            = 5 +1;  % Minimum Allowable UDREI (matlab runs 1:16 not 0:15 as in the MOPS)

UDRE_DW_ELEV =  7*pi/180;    % radians below which to start deweighting
UDRE_DW1     =  0.238060528495;
UDRE_DW2     = -1.17983911751;
UDRE_DW3     = -1.055141155463;

UDRE_MT28_MIN_CNMP2 = (0.45)^2;   % minimum CNMP value for MT28 in m

UDRE_PHI_MU     = 0.04219;           % Mean for min phi function
UDRE_PHI_SIG    = 0.4088067;         % Sigma for min phi function
UDRE_PHI_MIN_IN = 0.071364137584103; % Min input value to use above phi func.

UDRE_CLA_GPS = 1.1;  % critical limit adjustment
UDRE_K1 = 1.5;       % antenna bias parameters
UDRE_K2 = 0.03;
UDRE_K3 = 0.20;
UDRE_K4 = 21;
UDRE_K5 = 1.3;
UDRE_K6 = -0.015;
UDRE_K7 = 0.36;
UDRE_K8 = 47;










