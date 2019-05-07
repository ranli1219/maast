function [dCov,scalef,maxrat] = mt28cov(satxyz,satxyzdot,G,sig2_mon)

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

% Modified 2001 Nov01 By Todd Walter

global UDRE_SIG_SV_H UDRE_SIG_SV_C UDRE_SIG_SV_L UDRE_SIG_SV_T ...
        UDRE_SIG_DOT_H UDRE_SIG_DOT_C UDRE_SIG_DOT_L UDRE_T_PREDEG ...
        UDRE_MT28_MIN_CNMP2
% note: check signs of G matrix!

% a priori weighting matrix
inv_P0 = diag(1./[UDRE_SIG_SV_H^2 UDRE_SIG_SV_C^2 UDRE_SIG_SV_L^2 UDRE_SIG_SV_T^2]);
h = satxyz./norm(satxyz);
c = cross(satxyzdot./norm(satxyzdot),h);
l = cross(h,c);

H = [h;c;l];
H44 = eye(4);
H44(1:3,1:3) = H;

e=ones(size(sig2_mon));
sig2_mon=max(sig2_mon, e*UDRE_MT28_MIN_CNMP2);
W = diag(e./sig2_mon);
Cov = inv( H44'*inv_P0*H44 + G'*W*G );

% NORMALIZE Covariance with min projection
R = chol(Cov);
min_proj = R(4,4)^2;
nCov = Cov/min_proj;

% PREDEGRADE covariance
sig_hdot = UDRE_SIG_DOT_H*UDRE_T_PREDEG;
sig_cdot = UDRE_SIG_DOT_C*UDRE_T_PREDEG;
sig_ldot = UDRE_SIG_DOT_L*UDRE_T_PREDEG;
Pdelta = diag([sig_hdot sig_cdot sig_ldot].^2);
PEdelta = H'*Pdelta*H;


% calculate min projection los according to John Angus (02/23/01)
v0(3) = R(3,4)/R(3,3);
v0(2) = -(-R(2,4)*R(3,3)+R(2,3)*R(3,4))/(R(3,3)*R(2,2));
v0(1) = -(-R(1,4)*R(3,3)*R(2,2) + (R(2,4)*R(3,3)-R(2,3)*R(3,4))*R(1,2) + ...
          R(3,4)*R(2,2)*R(1,3))/(R(3,3)*R(2,2)*R(1,1));
PEv0 = PEdelta*v0';
nCov_star = [PEdelta, PEv0; ...
             PEv0' , v0*PEv0 ] / min_proj;
pdCov = nCov + nCov_star;

% DISCRETIZE covariance
[dCov,scalef] = disccov(pdCov);

maxrat = find_ja_maxratio(dCov,pdCov);




