function sig2 = af_geo_cnmpadd(del_t, ~)
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
%AF_GEO_CNMPADD calculate psueudorange confidence (variance) 
%               per GEO_CNMP ADD A014-008B
%SIG2=AF_WRSGEO_CNMPADD(DEL_T, EL)
%   DEL_T is the track time in seconds since last cycle slip
%   EL is the elevation angle, but is not used in this version
%   SIG2 is the psueudorange confidence (variance) in meters^2 
%   per GEO_CNMP ADD A014-008B
%
% SEE ALSO INIT_GEO_CNMP

%based on sigma_CNMP by Pete Schloss
%created 30 August, 2013 by Todd Walter

global GEO_CNMP_A0 GEO_CNMP_TL GEO_CNMP_TL2;
global GEO_CNMP_ME_TL GEO_CNMP_FLOOR GEO_CNMP_K GEO_CNMP_SIG_CARRIER;
global GEO_CNMP_A GEO_CNMP_B

mean_err=NaN(size(del_t));

idx=find(del_t < GEO_CNMP_TL/4);
if (~isempty(idx))
  mean_err(idx) = GEO_CNMP_ME_TL(del_t(idx));
end
idx=find((del_t>=GEO_CNMP_TL/4) & (del_t<GEO_CNMP_TL2));
if (~isempty(idx))
  mean_err(idx)=GEO_CNMP_A0./(2*pi*del_t(idx)/GEO_CNMP_TL);
end
idx=find(del_t>=GEO_CNMP_TL2);
if (~isempty(idx))
  mean_err(idx)=max(GEO_CNMP_FLOOR,(GEO_CNMP_A-(GEO_CNMP_B*(del_t(idx)-GEO_CNMP_TL2)))); 
end

sig2 = (mean_err/GEO_CNMP_K + GEO_CNMP_SIG_CARRIER).^2;

