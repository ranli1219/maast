function psi_max = mt28geo_cov2(G,sig2_mon, L)

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

% Created 2013 Aug 30 By Todd Walter

global GEO2_SIG_SV_EPH GEO2_SIG_SV_T GEO2_MT28_MIN_CNMP2

% note: check signs of G matrix!

% a priori weighting matrix
inv_P0 = diag(1./[GEO2_SIG_SV_EPH^2 GEO2_SIG_SV_EPH^2 GEO2_SIG_SV_EPH^2 GEO2_SIG_SV_T^2]);

e=ones(size(sig2_mon));
sig2_mon=max(sig2_mon, e*GEO2_MT28_MIN_CNMP2);
W = diag(e./sig2_mon);
P = inv(inv_P0 + G'*W*G);

% NORMALIZE Covariance with min projection
M = chol(P);
m44 = M(4,4)^2;
P_normalized = P/m44;

psi_max=find_psimax(P_normalized, L);






