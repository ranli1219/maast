function i = phmi_eq(sig2_cp, sig2_mon, mu_ant, alpha, lambda)

%*************************************************************************
%*     Copyright c 2004 The board of trustees of the Leland Stanford     *
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

global UDRE_KHMI UDRE_KFA UDRE_PHMI MOPS_UDRE UDRE_CLA_GPS   
global UDRE_SIG2_SV19 UDRE_SIG2_L1L2 UDRE_SIG2_OBAD   
global UDRE_PHI_MU UDRE_PHI_SIG UDRE_PHI_MIN_IN 
global MOPS_SIG_UDRE MOPS_UDREI_NM  % from init_mops

sig_obad = sqrt(UDRE_SIG2_OBAD);
sig_mon = sqrt(sig2_mon);

%initial value of udre_tbb
udre_tbb=sqrt(sig2_cp+UDRE_SIG2_SV19+UDRE_SIG2_L1L2)+sig_obad;
idx=find(MOPS_SIG_UDRE >= udre_tbb);
i=idx(1)-1;

%iterate over UDREIs
Pcalc = 1;
while Pcalc > UDRE_PHMI & i < MOPS_UDREI_NM
  i=i+1;
  udre_tbb=MOPS_SIG_UDRE(i);
  udre_star=sqrt((udre_tbb-sig_obad)^2-UDRE_SIG2_SV19-UDRE_SIG2_L1L2);

  K_a= (UDRE_KHMI*udre_star/lambda - UDRE_KFA*alpha.*sqrt(sig2_cp + sig2_mon) ...
                                   - UDRE_CLA_GPS - mu_ant)./sig_mon;

  %evaluate phi function
  pmu  = zeros(size(sig2_mon));
  psig =  ones(size(sig2_mon));
  below=find(K_a < UDRE_PHI_MIN_IN);
  if(~isempty(below))
    pmu(below)  = UDRE_PHI_MU;
    psig(below) = UDRE_PHI_SIG;
  end
  phi_a= 1-normcdf(K_a,pmu,psig);

  % N-1 constraint
  [phi_min,idx]=min(phi_a);
  phi_a(idx)=1;

  %calculate PI product
  Pcalc=prod(2*phi_a - phi_a.^2);
end

if(Pcalc > UDRE_PHMI)
  warning(' unable to meet Phmi requirement');
end



