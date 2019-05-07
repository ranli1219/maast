function init_geo2_osp()
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
%function init_geo_osp()
% 
% assigns values used by the UDRE ADD for Geostationary Satellites
%
% References: 	UDRE ADD Version 8

% Created 2003 Mar19 by Todd Walter
% Added Antenna Bias 2004 May21 by Todd Walter
% Updated to Release 3A1 2013 Aug 30 by Todd Walter

global GEO2_MT28_MIN_CNMP2 GEO2_SIG2_ROT100 GEO2_SIG2_WREL1L2 GEO2_SIG2_WRECLK
global GEO2_SIG2_SATL1L2 GEO2_PHMI GEO2_KHMI GEO2_KFA 
global GEO2_MIN_NWRS GEO2_UDREI_FLOOR
global GEO2_PHI_MIN_IN GEO2_PHI_MAX_IN 
global GEO2_K1 GEO2_K2 GEO2_K3 GEO2_K4 GEO2_CLA
global GEO2_DW1 GEO2_DW2 GEO2_DW3 GEO2_DW_ELEV GEO2_CEIL_EFF_GIVE
global GEO2_SIG_SV_EPH GEO2_SIG_SV_T GEO2_ALPHA_TROP

GEO2_MT28_MIN_CNMP2 = (.45)^2;
GEO2_SIG2_ROT100    = (0.8818)^2;        % updated R3
GEO2_SIG2_SATL1L2   = (.192/3.29)^2;
GEO2_SIG2_WREL1L2   = (.229/3.29)^2;    % updated R3
GEO2_SIG2_WRECLK    = (0.35/3.29)^2;

GEO2_SIG_SV_EPH  = 400000;
GEO2_SIG_SV_T    = 100;

GEO2_DW1     = 0.238060528495;  % First de-weighting OSP	
GEO2_DW2     = -1.17983911751;  % Second de-weighting OSP	
GEO2_DW3     = -1.055141155463; % Third de-weighting OSP
GEO2_DW_ELEV = 7*pi/180;        % Max elevation to apply deweighting

GEO2_ALPHA_TROP = (1/0.5)^2 - 1; %tropospheric adjustment for threshold

GEO2_CEIL_EFF_GIVE = 6/3.29;

GEO2_PHMI      = 4.5e-10;    % prob of HMI allocation to GEO2 monitor
GEO2_KHMI      = 5.33;       % matches MOPS KV_PA
GEO2_KFA       = 3.29/2;     % set empirically to trade avail. and continuity

%GEO2_UNOB_BIAS = 4.0;  %This term now defined in init_geo_cnmpadd

GEO2_MIN_NWRS          = 3;  % Minimum number of visible WRS

GEO2_UDREI_FLOOR       = 11;  % Minimum Allowable UDREI 7.5 m


GEO2_PHI_MIN_IN = 0.1;               % Min input value to use above phi func. % updated R3
GEO2_PHI_MAX_IN = 5;

GEO2_CLA = 0;      % critical limit adjustment

GEO2_K1 = 1.5;       % antenna bias parameters
GEO2_K2 = 0.03;
GEO2_K3 = 0.20;
GEO2_K4 = 21;








