function sig2_mvtec = sig2_max_vtec(t, delay, mag_lat, Kigd_fac, ceiling_P, ...
                                    ceiling_T, ceiling_lat)
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
%SIG2_MAX_VTEC returns the variance for the iono ceiling function
%
%  SIG2_MVTEC = SIG2_MAX_VTEC(T, DELAY, MAG_LAT)
%  Returns the variance corresponding to the ceiling funtion in the GIVE ADD
%  T is the local time of day at the IGP in seconds (between 0 and 86400)
%  DELAY is the delay estimate at the IGP in meters
%  MAG_LAT is the magnetic latitude of the IGP in degrees
%
%   See also : INIT_GIVE_OSP GIVE_ADD

%   TWalter 28 Mar 01


if mag_lat>ceiling_lat
  lat_idx=2;
else
  lat_idx=1;
end

if(t<=ceiling_T(lat_idx,1) | t>=ceiling_T(lat_idx,4))
  mvtec=ceiling_P(lat_idx,1);
elseif(t>ceiling_T(lat_idx,1) & t<ceiling_T(lat_idx,2))
  mvtec=ceiling_P(lat_idx,2);
elseif(t>ceiling_T(lat_idx,2) & t<ceiling_T(lat_idx,3))
  mvtec=ceiling_P(lat_idx,3);
else
  mvtec=ceiling_P(lat_idx,4);
end

sig2_mvtec=(Kigd_fac*max([delay mvtec-delay]))^2;

