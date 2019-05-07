function [sig2_igp, mu_igp, sig2_upm, r, rcm, idx, delay, chi2] ...
                           = igp_kriging(xyz_igp, igp_en_hat, ...
                                         xyz_ipp, mu_me, M, ...
                                         sig2_nom_decorr, ...
                                         sig2_tot_decorr, d_decorr, ...
                                         d_igp_corner, Ivpp)
%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
%IGP_PLANE returns the covariance matrix for a planar fit centered on the IGPs
%
%  [SIG2_IGP, MU_IGP, SIG2_UPM, R, RCM, IDX, DELAY, CHI2] = 
%                              IGP_PLANE(XYZ_IGP, IGP_EN_HAT, XYZ_IPP, M, IVPP)
%  XYZ_IGP are the XYZ ECEF coordinates of the IGP
%  IGP_EN_HAT are the unit vectors in the east and north directions (the east
%             vector is in the first 3 columns and the north vector is in 
%             columns 4 through 6
%  IGP_CORNERDEN are the corners for evaluating the projected antenna bias
%  XYZ_IPP are the XYZ ECEF coordinates of all the IPPs
%  MU_ME is maximum antenna bias on each line of sight in meters
%  M is the covariance matrix of the the IPPs used to form the fit
%  IVPP are the vertical delay for the IPPs (used for truth processing)
%  The return values are:
%  SIG2_IGP is the variance for the IGP delay fit
%  MU_IGP is the maximum antenna bias effect on the fit
%  SIG2_UPM is teh UPM variance from the fit, used for GEO UDRE
%  R is the radius used for the fit
%  RCM is the relative centroid (weighted average of IPP positions / R)
%  IDX are the indices of the IPPs used in the fit
%  DELAY is the vertical delay estimate for the IGP (truth only)
%  CHI2 is the chi-squared value for the fit (truth only)
%
%   See also: GIVE_ADD1 INIT_GIVEADD1_OSP

%   created 18 May 2004 by Todd Walter
%   updated 29 June 2006 by Todd Walter: added antenna bias term
%   updated 19 Nov. 2009 by Todd Walter: removed rcmx
%   updated 23 August 2013 by Todd Walter to align with Release 3 ADD

global MOPS_NOT_MONITORED
global GIVEK_N_MIN GIVEK_N_PTS GIVEK_RMAX GIVEK_RMIN  ...
       GIVEK_CHI2_LOWER GIVEK_CHI2_LWR_UPM GIVEK_CHI2_NOM GIVEK_RIRREG2_FLOOR
global TRUTH_FLAG

%initialize values
sig2_igp = MOPS_NOT_MONITORED;
mu_igp = 0;
sig2_upm = MOPS_NOT_MONITORED;
rcm=1;
delay=NaN;
chi2=NaN;

n_ipps = size(xyz_ipp,1);

%calculate distance from IGP to each IPP
e=ones(n_ipps,1);
del_xyz=xyz_ipp-e*xyz_igp;
dist2=sum(del_xyz'.^2)';

% find at least 30 IPPs between max and min radii
%  try max radius first
r=GIVEK_RMAX;
idx=find(dist2 <= GIVEK_RMAX^2);
n_fit=length(idx);
%  do we meet the minimum requirement?
if (n_fit < GIVEK_N_MIN)
%  warning(['Only ', num2str(n_fit), ' IPPs fell within radius ',...
%          num2str(GIVE_RMAX/1e3), ' km of IGP']);
  return
end
%  can we use a smaller radius?
if (n_fit > GIVEK_N_PTS)
  min_idx=find(dist2(idx) <= GIVEK_RMIN^2);
  n_fit=length(min_idx);
%    check minimum radius
  if(n_fit >= GIVEK_N_PTS)
    r=GIVEK_RMIN;
    idx=idx(min_idx);
  else
    [temp, fit_idx]=sort(dist2(idx));
    r=sqrt(temp(GIVEK_N_PTS));
    idx=idx(fit_idx(1:GIVEK_N_PTS));
    n_fit=GIVEK_N_PTS;
  end
end



%build the G, C, and W matrices
V = del_xyz(idx,:)*1e-6;
             
G=[e(idx) sum((V.*(e(idx)*igp_en_hat(:,1:3))),2)...
                 sum((V.*(e(idx)*igp_en_hat(:,4:6))),2)];

%calculate distance from each IPP to all other IPPs
D = reshape(sqrt(sum((repmat(V,n_fit,1) - ...
      reshape(repmat(V(:)',n_fit,1),n_fit^2,3)).^2,2)),n_fit,n_fit);             
d = sqrt(sum(V.^2,2));

C = (sig2_tot_decorr - sig2_nom_decorr)*exp(-D/d_decorr);
C_min = min(C)';
C = C + (sig2_tot_decorr - C(1,1))*eye(n_fit);

c = (sig2_tot_decorr - sig2_nom_decorr)*exp(-d/d_decorr);

delta = (sig2_tot_decorr - sig2_nom_decorr)*(1-exp(-d_igp_corner/d_decorr)); 

W=inv(M(idx,idx) + C);
  
%calculate covariance matrix
cov0 = inv(G'*W*G); 

P = G*(cov0)*G'*W;

w = (W - W*P)*c + W*G*(cov0)*([1 0 0]');

if TRUTH_FLAG
    %calculate vertical delay value
    delay = w'*Ivpp(idx);

    %calculate chi-square value of noiseless supertruth
    chi2_raw=Ivpp(idx)'*W*(eye(n_fit) - P)*Ivpp(idx);
    Rnoise = 1;
else
    
    Rnoise = prod(1+diag(M(idx,idx))./(diag(C) - C_min))^(1/n_fit);
    chi2_raw = 0.65*GIVEK_CHI2_NOM(n_fit-3)/Rnoise;
end

chi2 = chi2_raw*Rnoise;

Rirreg2=chi2/GIVEK_CHI2_LOWER(n_fit-3);
Rirreg2=max(GIVEK_RIRREG2_FLOOR,Rirreg2);

sig2_igp = Rirreg2*(w'*C*w -2*w'*c + sig2_tot_decorr + delta) + ...
                         w'*M(idx,idx)*w;

%calculate antenna bias term
mu_igp = sum(abs(w).*mu_me(idx));

%calculate UPM GIVE value for GEO UDRE
Rirreg_realistic2 = chi2/GIVEK_CHI2_LWR_UPM(n_fit-3);
Rirreg_realistic2 = max(GIVEK_RIRREG2_FLOOR,Rirreg_realistic2);
sig2_upm = Rirreg_realistic2*(w'*C*w -2*w'*c + sig2_tot_decorr) + ...
                         w'*M(idx,idx)*w;

%calculate Relative Centroid Metric
We=e(idx)./(diag(W));

rcv=sum((del_xyz(idx,:).*[We We We]))'/sum(We);
rcm=sqrt(sum(rcv.^2))/r;


       
       
