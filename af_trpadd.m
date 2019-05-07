function sig2_trop = af_wrstrpadd(El)
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
% tropo confidence bound (variance) from UDRE ADD
% El can be any n-by-m matrix of elevations.
% Returns tropo variance for each elevation element.  
% NaN is returned for NaN elevations.

global TROP_SIG_VERT
 
idxvis = find(~isnan(El));
sig2_trop = repmat(NaN,size(El));
% from ADD
sig2_trop(idxvis) = (TROP_SIG_VERT./sin(El(idxvis))).^2;                


