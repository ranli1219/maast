function sig2_undersamp=undersamp_threat(radius, rel_centroid, critical_pts)

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
%UNDERSAMP_THREAT calculates the undersampled threat variance per the ADD
%
%  SIG2_UNDERSAMP=UNDERSAMP_THREAT(RADIUS, REL_CENTROID)
%  RADIUS (km) is the radius encompassing the IPPs used in the planar fit
%  REL_CENTROID is the weighted centroid of the IPP positions used in the fit
%  divided by the radius.  it is dimensionless.
%
% SEE ALSO: INIT_GIVE_OSP GIVE_ADD

%created 28 Mar 2001 by Todd Walter
%modified 05 Nov 2001 by Todd Walter


idx=find(radius >= critical_pts(:,2) & ...
         rel_centroid >= critical_pts(:,1));

if(~isempty(idx))
  sig2_undersamp=max(critical_pts(idx,3)).^2;
else
  sig2_undersamp=max(critical_pts(:,3)).^2;
end



