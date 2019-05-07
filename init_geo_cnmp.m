function init_geo_cnmp()
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
%SIG2_GEO_CNMP calculate psueudorange confidence (variance) per GEO_CNMP ADD A014-008B
% initializes GEO_CNMP OSPs
%   per CNMP ADD A014-008B
% SEE ALSO: SIG2_GEO_CNMP

%based on sigma_CNMP by Pete Schloss
%created 30 August, 2013 by Todd Walter


global GEO_CNMP_A0 GEO_CNMP_TL GEO_CNMP_TL2;
global GEO_CNMP_ME_TL GEO_CNMP_FLOOR GEO_CNMP_K GEO_CNMP_SIG_CARRIER;
global GEO_CNMP_A GEO_CNMP_B GEO_CNMP_UNOB_BIAS_INIT
global GEO_CNMP_UNOB_BIAS GEO_CNMP_UNOB_BIAS_INIT_CNMP_VAL

GEO_CNMP_A0          = 30;
GEO_CNMP_TL          = 86400;
GEO_CNMP_TL1         = GEO_CNMP_TL/4;
GEO_CNMP_TL2         = 43200;
GEO_CNMP_TL3         = 259200;
GEO_CNMP_FLOOR       = 1.0;
GEO_CNMP_K           = 3.29;
GEO_CNMP_SIG_CARRIER = 0.03;

GEO_CNMP_A           = GEO_CNMP_A0/(2*pi*GEO_CNMP_TL2/GEO_CNMP_TL);
GEO_CNMP_B           = (GEO_CNMP_A-GEO_CNMP_FLOOR)/(GEO_CNMP_TL3-GEO_CNMP_TL2);

t=1:GEO_CNMP_TL1;
GEO_CNMP_ME_TL       = GEO_CNMP_A0*sin(2*pi*t/GEO_CNMP_TL)./(2*pi*t/GEO_CNMP_TL);

GEO_CNMP_UNOB_BIAS_INIT = 4.0;
GEO_CNMP_UNOB_BIAS      = 3.5;
GEO_CNMP_UNOB_BIAS_INIT_CNMP_VAL = af_geo_cnmpadd(129600);
