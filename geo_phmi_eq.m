function i = geo_phmi_eq(sig2_cp, sig2_mon, mu_ant, lambda2_wrs)

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
% function i = geo_phmi_eq(sig2_cp,sig2_mon)
% Calculates udre index based on phmi equation.
% Output:   i           - one-based udre index 
%
% Inputs:   sig_cp [m]  - udre from corrections processor
%           sig_mon [m] - monitor noise sigma 

% Modified 2003 Mar20 By Todd Walter
% Added Antenna Bias 2004 May21 by Todd Walter

global GEO_KHMI GEO_KFA GEO_PHMI GEO_SIG2_SATL1L2 GEO_SIG2_GIVEFLOOR
global GEO_PHI_MU GEO_PHI_SIG GEO_PHI_MIN_IN GEO_PHI_MAX_IN GEO_UNOB_BIAS
global GEO_CLA
global MOPS_SIG2_UDRE MOPS_UDREI_NM  % from init_mops

sig_mon=sqrt(sig2_mon);

%initial value of udre_tbb
udre_tbb=sig2_cp+GEO_SIG2_SATL1L2;
idx=find(MOPS_SIG2_UDRE >= udre_tbb);
i=idx(1)-1;

%iterate over UDREIs
Pcalc = 1;
while Pcalc > GEO_PHMI & i < MOPS_UDREI_NM
  i=i+1;
  udre_tbb=MOPS_SIG2_UDRE(i);
  sig2_udre_star=udre_tbb-GEO_SIG2_SATL1L2;

  K_a = (GEO_KHMI*sqrt(lambda2_wrs*sig2_udre_star + GEO_SIG2_GIVEFLOOR) - ...
           GEO_KFA*sqrt(sig2_cp + sig2_mon) - GEO_UNOB_BIAS - GEO_CLA - ...
           mu_ant)./sig_mon;
  
  % max on the phi function
  above=find(K_a > GEO_PHI_MAX_IN);
  if(~isempty(above))
    K_a(above)  = GEO_PHI_MAX_IN;
  end
  %evaluate phi function
  pmu  = zeros(size(sig2_mon));
  psig =  ones(size(sig2_mon));
  below=find(K_a < GEO_PHI_MIN_IN);
  if(~isempty(below))
    pmu(below)  = GEO_PHI_MU;
    psig(below) = GEO_PHI_SIG;
  end
  phi_a= 1-normcdf(K_a,pmu,psig);

  % N-1 constraint
  [phi_min,idx]=min(phi_a);
  phi_a(idx)=1;

  %calculate PI product
  Pcalc=prod(2*phi_a - phi_a.^2);
end

if(Pcalc > GEO_PHMI)
  warning(' unable to meet Phmi requirement');
end



