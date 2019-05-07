function init_geo_osp()
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
%function init_geo_osp()
% 
% assigns values used by the UDRE ADD for Geostationary Satellites
%
% References: 	UDRE ADD Version 8

% Created 2003 Mar19 by Todd Walter
% Added Antenna Bias 2004 May21 by Todd Walter

global GEO_SIG2_CNMPFLOOR GEO_SIG2_ROT100 GEO_SIG2_WREL1L2 GEO_SIG2_WRECLK
global GEO_SIG2_SATL1L2 GEO_PHMI GEO_KHMI GEO_KFA GEO_SIG2_GIVEFLOOR
global GEO_MIN_NWRS GEO_MIN_NWRSCONUS GEO_UDREI_FLOOR_NOMINCONUS GEO_UDREI_FLOOR
global GEO_PHI_MU GEO_PHI_SIG GEO_PHI_MIN_IN GEO_PHI_MAX_IN GEO_UNOB_BIAS
global GEO_K1 GEO_K2 GEO_K3 GEO_K4 GEO_CLA
		
		
GEO_SIG2_CNMPFLOOR = (.334)^2;
GEO_SIG2_ROT100 = (0.3508)^2;
GEO_SIG2_SATL1L2 = (.192)^2;
GEO_SIG2_WREL1L2 = (.006)^2;
GEO_SIG2_WRECLK = (0.35)^2;


GEO_SIG2_GIVEFLOOR = (1.2*3.0/3.29)^2;


GEO_PHMI      = 4.5e-10;    % prob of HMI allocation to GEO monitor
GEO_KHMI      = 5.33;       % matches MOPS KV_PA
GEO_KFA       = 3.29/2;     % set empirically to trade avail. and continuity

GEO_UNOB_BIAS = 5.0;

GEO_MIN_NWRS          = 3;  % Minimum number of visible WRS
GEO_MIN_NWRSCONUS     = 1;  % Minimum number of visible CONUS WRS
GEO_UDREI_FLOOR_NOMINCONUS = 7; % Udrei floor for enough wrs but not 
                                 %  enough conus wrs
GEO_UDREI_FLOOR       = 7;  % Minimum Allowable UDREI

GEO_PHI_MU     = 0.04219;           % Mean for min phi function
GEO_PHI_SIG    = 0.4088067;         % Sigma for min phi function
GEO_PHI_MIN_IN = 0.071364137584103; % Min input value to use above phi func.
GEO_PHI_MAX_IN = 5;

GEO_CLA = 0;      % critical limit adjustment
GEO_K1 = 1.5;       % antenna bias parameters
GEO_K2 = 0.03;
GEO_K3 = 0.20;
GEO_K4 = 21;









