function i = geo_phmi_eq2(sig2_cp, sig2_mon, sig2_trop, mu_ant, ...
                          alpha_dw, psi_max, min_uire2_sp, geo_unob_bias)

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
% function i = geo_phmi_eq(sig2_cp,sig2_mon)
% Calculates udre index based on phmi equation.
% Output:   i           - one-based udre index 
%
% Inputs:   sig_cp [m]  - udre from corrections processor
%           sig_mon [m] - monitor noise sigma 

% Modified 2003 Mar20 By Todd Walter
% Added Antenna Bias 2004 May21 by Todd Walter
% Updated to Release 8/9 Nov 20, 2009 by Todd Walter
% Updated to R3A Aug 30, 2013  by Todd Walter

global GEO2_KHMI GEO2_KFA GEO2_PHMI GEO2_SIG2_SATL1L2 
global GEO2_PHI_MIN_IN GEO2_PHI_MAX_IN
global GEO2_CLA GEO2_ALPHA_TROP
global MOPS_SIG2_UDRE MOPS_UDREI_NM  % from init_mops

sig_mon=sqrt(sig2_mon);

%initial value of udre_tbb
udre_tbb=sig2_cp+GEO2_SIG2_SATL1L2;
idx=find(MOPS_SIG2_UDRE >= udre_tbb);
i=idx(1)-1;

%iterate over UDREIs
Pcalc = 1;
while Pcalc > GEO2_PHMI && i < MOPS_UDREI_NM
  i=i+1;
  udre_tbb=MOPS_SIG2_UDRE(i);
  sig2_udre_star=udre_tbb-GEO2_SIG2_SATL1L2;

  K_a = (GEO2_KHMI*sqrt(sig2_udre_star/psi_max + min_uire2_sp) - ...
           GEO2_KFA*alpha_dw.*sqrt(sig2_cp + sig2_mon + GEO2_ALPHA_TROP*sig2_trop)...
           - geo_unob_bias - GEO2_CLA - mu_ant)./sig_mon;
  
  % max on the phi function
  above=find(K_a > GEO2_PHI_MAX_IN);
  if(~isempty(above))
    K_a(above)  = GEO2_PHI_MAX_IN;
  end
  %evaluate phi function
  phi_a= 1-normcdf(K_a);
  phi_a(K_a < GEO2_PHI_MIN_IN) = 1;

  % N-1 constraint
  [~,idx]=min(phi_a);
  phi_a(idx)=1;

  %calculate PI product
  Pcalc=prod(2*phi_a - phi_a.^2);
end

if(Pcalc > GEO2_PHMI)
  warning(' unable to meet Phmi requirement');
end



