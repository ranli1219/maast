function sig2=sig2_cnmp(del_t,el)
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
%SIG2=SIG2_CNMP(DEL_T, EL)
%   DEL_T is the track time in seconds since last cycle slip
%   EL is the elevation angle, but is not used in this version
%   SIG2 is the psueudorange confidence (variance) in meters^2 
%   per CNMP ADD A014-008B
%
% SEE ALSO INIT_CNMP

%based on sigma_CNMP by Pete Schloss
%created 20 April, 2001 by Todd Walter

global CNMP_A0 CNMP_TL CNMP_TL1 CNMP_TL2 CNMP_TL3;
global CNMP_ME_TL CNMP_FLOOR CNMP_K CNMP_SIG_CARRIER;
global CNMP_A CNMP_B

mean_err=repmat(NaN,size(del_t));

idx=find(del_t < CNMP_TL/4);
if (~isempty(idx))
  mean_err(idx) = CNMP_ME_TL(del_t(idx));
end
idx=find((del_t>=CNMP_TL/4) & (del_t<CNMP_TL2));
if (~isempty(idx))
  mean_err(idx)=CNMP_A0./(2*pi*del_t(idx)/CNMP_TL);
end
idx=find(del_t>=CNMP_TL2);
if (~isempty(idx))
  mean_err(idx)=max(CNMP_FLOOR,(CNMP_A-(CNMP_B*(del_t(idx)-CNMP_TL2)))); 
end

sig2 = (mean_err/CNMP_K + CNMP_SIG_CARRIER).^2;

   
