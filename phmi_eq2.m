function i = phmi_eq2(sig2_cp, sig2_mon, sig2_trop, mu_ant, alpha_dw, lambda)

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
% function i = phmi_eq(sig2_cp,sig2_mon)
% Calculates udre index based on phmi equation.
% Output:   i           - one-based udre index 
%
% Inputs:   sig_cp [m]  - udre from corrections processor
%           sig_mon [m] - monitor noise sigma 

% Modified 2001 Nov01 By Todd Walter
% Added Antenna Bias Parameters 2004 May21 by Todd Walter
% Updated to Release 8/9 Nov 20, 2009 by Todd Walter
% Updated to Release 3A1 Sept. 4, 2013 by Todd Walter

global UDRE2_KHMI UDRE2_KFA UDRE2_PHMI UDRE2_CLA_GPS   
global UDRE2_SIG2_L1L2 UDRE2_ALPHA_TROP
global UDRE2_PHI_LOWER_BND UDRE2_PHI_CUTOFF
global MOPS_UDREI_NM MOPS_SIG_UDRE % from init_mops

sig_mon = sqrt(sig2_mon);

%initial value of udre_tbb
udre_tbb=sqrt(sig2_cp+UDRE2_SIG2_L1L2);
idx=find(MOPS_SIG_UDRE >= udre_tbb);
i=idx(1)-1;

%iterate over UDREIs
Pcalc = 1;
while Pcalc > UDRE2_PHMI && i < MOPS_UDREI_NM
  i=i+1;
  udre_tbb=MOPS_SIG_UDRE(i);
  udre_star=sqrt(udre_tbb^2-UDRE2_SIG2_L1L2);

  K_a= (UDRE2_KHMI*udre_star/lambda - UDRE2_KFA*alpha_dw.*...
               sqrt(sig2_cp + sig2_mon + UDRE2_ALPHA_TROP*sig2_trop) ...
                                   - UDRE2_CLA_GPS - mu_ant)./sig_mon;

  %evaluate phi function
  phi_a= 1-normcdf(K_a);
  phi_a(K_a < UDRE2_PHI_LOWER_BND) = 1;
  phi_a(K_a > UDRE2_PHI_CUTOFF) = 1-normcdf(UDRE2_PHI_CUTOFF);

  % N-1 constraint
  [~,idx]=min(phi_a);
  phi_a(idx)=1;

  %calculate PI product
  Pcalc=prod(2*phi_a - phi_a.^2);
end

if(Pcalc > UDRE2_PHMI)
  warning(' unable to meet Phmi requirement');
end



