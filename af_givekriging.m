function igpdata = af_givekriging(~, igpdata, wrsdata, satdata, ...
                                  wrs2satdata, ~)
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
%GIVE_ADD7 returns the give variances for each IGP
%
%  GIVEI = GIVE_ADD(T, XYZ_IGP, IGP_EN_HAT, IGP_CORNER_DEN, IGP_LON, 
%                   IGP_MAG_LAT, XYZ_IPP, SIG2_IPP, WRS_IDX)
%  T is the GPS time of day in seconds
%  XYZ_IGP (n_igp,3) are the XYZ ECEF coordinates of all the IGPs (m)
%  IGP_EN_HAT (n_igp,6) are the unit vectors in the east and north direction
%             (the east vector is in the first 3 columns and the north vector
%             is in columns 4 through 6)
%  IGP_CORNER_DEN (4,3,n_igp) are the [1 delta_east(Mm) delta_north(Mm)]
%                 vectors for the cell centers surrounding the IGP
%  IGP_LON (n_igp,1) is the IGP longitude in degrees
%  IGP_MAG_LAT (n_igp,1) are the magnetic latitudes for the IGPs (degrees)
%  XYZ_IPP (n_ipp,3) are the XYZ ECEF coordinates of all the IPPs (m)
%  SIG2_IPP (n_ipp,1) are the variances of the the IPPs used to weight the fit
%           (m^2)
%  WRS_IDX (n_ipp,1) are the WRS identifying numbers for the IPPs
%  The return values are:
%  GIVEI the discretized confidence variance for each IGP
%
%   See also: INIT_GIVEADD1_OSP IGP_PLANE1 UNDERSAMP_THREAT CHECKFOR2

%   Created 18 May 2004 by Todd Walter
%   updated 3 September 2013 by Todd Walter to align with Release 3 ADD


%get constant parameters
global CONST_F1 CONST_F2
global MOPS_SIG2_GIVE MOPS_GIVEI_NM
global  GIVEK_N_MIN ...
        GIVEK_KA_RATIO GIVEK_SIG2_L1L2REC GIVEK_SIG2_L1L2SAT ...
        GIVEK_GIVEI_FLOOR GIVEK_UNDERSAMP_PTS GIVEK_CHI2_THRESH ...
        GIVEK_MEX_SLOPE1 GIVEK_MEX_SLOPE2 GIVEK_MEX_LLAT1 GIVEK_MEX_LLAT2 ...
        GIVEK_ANT_A1 GIVEK_ANT_B1 GIVEK_ANT_C1 GIVEK_ANT_ALPHA1 ...
        GIVEK_ANT_A2 GIVEK_ANT_B2 GIVEK_ANT_C2 GIVEK_ANT_ALPHA2 ...
        GIVEK_UPMUNSAMP_PTS GIVEK_KHMI_KA
global GIVEK_SIG2_NOM GIVEK_SIG2_TOT GIVEK_D_DECORR ...
       GIVEK_D_IGP_CORNER5X5 GIVEK_D_IGP_CORNER5X10
    
global   COL_IGP_LL COL_IGP_XYZ COL_IGP_FLOORI ...
        COL_IGP_EHAT COL_IGP_NHAT  ...
        COL_IGP_GIVEI COL_IGP_MINMON COL_IGP_DELAY  ...
        COL_IGP_UPMGIVEI ...
        COL_USR_UID COL_USR_LL ...
        COL_U2S_UID COL_U2S_PRN   ...
          COL_U2S_EL  ...
        COL_U2S_SIG2L1MP COL_U2S_SIG2L2MP COL_U2S_IPPLL COL_U2S_IPPXYZ ...
         COL_U2S_IVPP   
global TRUTH_FLAG  TRIP_COUNT;

n_igp = size(igpdata,1);
xyz_igp = igpdata(:,COL_IGP_XYZ);
igp_en_hat = igpdata(:,[COL_IGP_EHAT,COL_IGP_NHAT]);
igp_lat = igpdata(:,COL_IGP_LL(1));

%find all IPPs above the elevation mask and above the cutoff line for
%mexico stations
keep_ipps = ones(size(wrs2satdata,1),1);
mex_wrs = find(wrsdata(:,COL_USR_LL(1))< 24.1);
if ~isempty(mex_wrs)
    mex_wrs = wrsdata(mex_wrs,COL_USR_UID);
    low_mex_ipps = find(ismember(wrs2satdata(:,COL_U2S_UID),mex_wrs) & ...
                        (wrs2satdata(:,COL_U2S_IPPLL(1)) < ...
                        GIVEK_MEX_SLOPE1*wrs2satdata(:,COL_U2S_IPPLL(2)) + ...
                        GIVEK_MEX_LLAT1 | wrs2satdata(:,COL_U2S_IPPLL(1)) < ...
                        GIVEK_MEX_SLOPE2*wrs2satdata(:,COL_U2S_IPPLL(2)) + ...
                        GIVEK_MEX_LLAT2));
    if ~isempty(low_mex_ipps)
        keep_ipps(low_mex_ipps) = 0;
    end
end
abv_mask = find(~isnan(wrs2satdata(:,COL_U2S_EL)) & keep_ipps);
xyz_ipp = wrs2satdata(abv_mask,COL_U2S_IPPXYZ);
nwrs = size(wrsdata,1);
nlos = size(wrs2satdata,1);
nsat = size(satdata,1);
wrs_idx = reshape(repmat(1:nwrs,nsat,1),nlos,1);
wrs_idx = wrs_idx(abv_mask);
inv_obl2 = ones(length(abv_mask),1)./...
              obliquity2(wrs2satdata(abv_mask,COL_U2S_EL));
            
            
M = diag((((CONST_F2^2/(CONST_F1^2-CONST_F2^2)).^2)*...
         ( wrs2satdata(abv_mask,COL_U2S_SIG2L1MP)+...
           wrs2satdata(abv_mask,COL_U2S_SIG2L2MP) )) ...
           .*inv_obl2 );
%fill in the off-diagonal elements
obl_mat=sqrt(inv_obl2*inv_obl2');
sat_idx=wrs2satdata(abv_mask,COL_U2S_PRN);
wrs_s=unique(wrs_idx);
for idx = 1:length(wrs_s);
    elem = find(wrs_idx==wrs_s(idx));
    M(elem,elem) = M(elem,elem) + GIVEK_SIG2_L1L2REC*obl_mat(elem,elem);
end
prn_s=unique(sat_idx);
for idx = 1:length(prn_s);
    elem = find(sat_idx==prn_s(idx));
    M(elem,elem) = M(elem,elem) + GIVEK_SIG2_L1L2SAT*obl_mat(elem,elem);
end

Ivpp = wrs2satdata(abv_mask,COL_U2S_IVPP);

%antenna bias terms
mu_me = (CONST_F2^2/(CONST_F1^2-CONST_F2^2))*sqrt(inv_obl2).*...
        (GIVEK_ANT_C1*(GIVEK_ANT_A1 + GIVEK_ANT_B1*...
         exp(-wrs2satdata(abv_mask,COL_U2S_EL)/GIVEK_ANT_ALPHA1)) + ...
         GIVEK_ANT_C2*(GIVEK_ANT_A2 + GIVEK_ANT_B2*...
         exp(-wrs2satdata(abv_mask,COL_U2S_EL)/GIVEK_ANT_ALPHA2)));

%initialize return values
igpdata(:,COL_IGP_MINMON)=zeros(n_igp,1);
givei=repmat(MOPS_GIVEI_NM,n_igp,1);
upm_givei=givei-2;

%loop over each Ionospheric Grid Point
for igp=1:n_igp
  %calculate the covariance matrix, radius, relative centroid metric and get
  %the indices for the IPPs used in the planar fit
  d_igp_corner = GIVEK_D_IGP_CORNER5X5;
  if igp_lat > 55
        d_igp_corner = GIVEK_D_IGP_CORNER5X10;
  end
  [sig2_igp, mu_igp, sig2_upm, r, rcm, idx, delay, chi2] = ...
                igp_kriging(xyz_igp(igp,:), igp_en_hat(igp,:), ...
                           xyz_ipp, mu_me, M, GIVEK_SIG2_NOM, ...
                           GIVEK_SIG2_TOT, GIVEK_D_DECORR, ...
                           d_igp_corner, Ivpp);
           
  if TRUTH_FLAG
      if delay < 0
        delay = 0;
      end
      igpdata(igp,COL_IGP_DELAY) = delay;
  end
                                                                     
  n_ipp=length(idx);
  %make sure we have the minimum number of points

  if(n_ipp>=GIVEK_N_MIN)
	%check for at least 2 WRSs with 2 IPPs
    if(checkfor2(wrs_idx(idx)))
	  % flag that this meets minimum monitoring requirements
	  igpdata(igp,COL_IGP_MINMON)=1;

      %calculate the undersampled threat
      sig2_threat = undersamp_threat(r/1000, rcm, GIVEK_UNDERSAMP_PTS);
      sig2_upm_threat = undersamp_threat(r/1000, rcm, GIVEK_UPMUNSAMP_PTS);
      
      %assemble GIVE 
      sig2_give=GIVEK_KA_RATIO^2*((sqrt(sig2_igp) + mu_igp/GIVEK_KHMI_KA)^2 + ...
                                     sig2_threat);

      %discretize the GIVE value
      if(sig2_give <= MOPS_SIG2_GIVE(MOPS_GIVEI_NM-1))
        temp=find(floor(MOPS_SIG2_GIVE/sig2_give));
        givei(igp)=temp(1);
        % check per IGP floor
        if givei(igp) < igpdata(igp, COL_IGP_FLOORI)
            givei(igp) = igpdata(igp, COL_IGP_FLOORI);
        end
        % UPM GIVE
        temp=find(floor(MOPS_SIG2_GIVE/(sig2_upm + sig2_upm_threat)));
        upm_givei(igp)=temp(1);
        
        %storm detector
        if chi2 >= GIVEK_CHI2_THRESH(n_ipp-3)
            TRIP_COUNT = TRIP_COUNT + 1;
            givei(igp) = MOPS_GIVEI_NM-1;
            upm_givei(igp) = MOPS_GIVEI_NM-1;
        end
      end   
    end
  end
end
%check global floor
below_floor=find(givei < GIVEK_GIVEI_FLOOR);
if(~isempty(below_floor))
  givei(below_floor)=GIVEK_GIVEI_FLOOR;
end

igpdata(:,COL_IGP_GIVEI) = givei;
igpdata(:,COL_IGP_UPMGIVEI) = upm_givei;
