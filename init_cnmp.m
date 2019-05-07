function init_cnmp()
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
%SIG2_CNMP calculate psueudorange confidence (variance) per CNMP ADD A014-008B
% initializes CNMP OSPs
%   per CNMP ADD A014-008B
% SEE ALSO: SIG2_CNMP

%based on sigma_CNMP by Pete Schloss
%created 20 April, 2001 by Todd Walter
% modified 31 April 2001 by Wyant Chan 
%       - added step intervals for time track

global CNMP_A0 CNMP_TL CNMP_TL1 CNMP_TL2 CNMP_TL3;
global CNMP_ME_TL CNMP_FLOOR CNMP_K CNMP_SIG_CARRIER;
global CNMP_A CNMP_B

CNMP_A0          = 10;
CNMP_TL          = 1600;
CNMP_TL1         = CNMP_TL/4;
CNMP_TL2         = 1525;
CNMP_TL3         = 12000;
CNMP_FLOOR       = 0.4;
CNMP_K           = 3.29;
CNMP_SIG_CARRIER = 0.03;

CNMP_A           = CNMP_A0/(2*pi*CNMP_TL2/CNMP_TL);
CNMP_B           = (CNMP_A-CNMP_FLOOR)/(CNMP_TL3-CNMP_TL2);

t=1:CNMP_TL1;
CNMP_ME_TL       = CNMP_A0*sin(2*pi*t/CNMP_TL)./(2*pi*t/CNMP_TL);

