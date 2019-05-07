function igpdata = af_giveadd(t, igpdata, wrsdata, satdata, wrs2satdata, ...
                              truth_data)
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
%GIVE_ADD returns the give variancs for each IGP
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
%   See also: INIT_GIVE_OSP IGP_PLANE SIG2_MAX_VTEC UNDERSAMP_THREAT CHECKFOR2

%   Created 28 Mar 2001 by Todd Walter
%   Modifications:
%       09 Apr 2001 Wyant Chan
%       1) Made wrapper around igp data to conform with IGPDATA format.
%       2) Return entire igpdata matrix instead of just the givei's.
%       3) Made wrapper to conform with WRS2SATDATA format
%       4) ipp data matrices runs through all los's, including the non
%           visible los's.  NaN's are used.
%modified 05 Nov 2001 by Todd Walter
%modified 16 Jun 2003 by Matt DeLand to implement Larry Spark's metric, IPP
%spread

%get constant parameters
global CONST_SEC_PER_DAY CONST_F1 CONST_F2
global MOPS_SIG2_GIVE MOPS_GIVEI_NM
global GIVE_SIG2_DECORR GIVE_R2_IRREG GIVE_N_MIN ...
        GIVE_SIG2_ROT0 GIVE_KHMI_KA GIVE_SIG2_L1L2REC GIVE_SIG2_L1L2SAT ...
        GIVE_GIVEI_FLOOR GIVE_UNDERSAMP_PTS GIVE_KIGD_FAC GIVE_CEILING_P...
        GIVE_CEILING_T GIVE_CEILING_LAT GIVE_CHI2_THRESH  GIVE_R2_DYN_IRREG
global COL_IGP_BAND COL_IGP_ID COL_IGP_LL COL_IGP_XYZ COL_IGP_ISCONUS ...
        COL_IGP_EHAT COL_IGP_NHAT COL_IGP_MAGLAT COL_IGP_CORNERDEN ...
        COL_IGP_GIVEI COL_IGP_MINMON COL_IGP_DELAY COL_IGP_MAX ...
        COL_IGP_UPMGIVEI  COL_IGP_BETA ...
        COL_U2S_UID COL_U2S_PRN COL_U2S_LOSXYZ COL_U2S_GXYZB ...
        COL_U2S_LOSENU COL_U2S_GENUB COL_U2S_EL COL_U2S_AZ COL_U2S_SIG2TRP ...
        COL_U2S_SIG2L1MP COL_U2S_SIG2L2MP COL_U2S_IPPLL COL_U2S_IPPXYZ ...
        COL_U2S_TTRACK0 COL_U2S_IVPP COL_U2S_MAX COL_U2S_INITNAN 
global TRUTH_FLAG RTR_FLAG RTR_MULT IPP_SPREAD_FLAG TRIP_COUNT;
      
n_igp = size(igpdata,1);
xyz_igp = igpdata(:,COL_IGP_XYZ);
igp_en_hat = igpdata(:,[COL_IGP_EHAT,COL_IGP_NHAT]);
igp_corner_den = reshape(igpdata(:,COL_IGP_CORNERDEN)',4,3,n_igp);
igp_lon = igpdata(:,COL_IGP_LL(2));
igp_mag_lat = igpdata(:,COL_IGP_MAGLAT);
abv_mask = find(~isnan(wrs2satdata(:,COL_U2S_EL)));
xyz_ipp = wrs2satdata(abv_mask,COL_U2S_IPPXYZ);
nwrs = size(wrsdata,1);
nlos = size(wrs2satdata,1);
nsat = size(satdata,1);
wrs_idx = reshape(repmat(1:nwrs,nsat,1),nlos,1);
wrs_idx = wrs_idx(abv_mask);
inv_obl2 = ones(length(abv_mask),1)./...
              obliquity2(wrs2satdata(abv_mask,COL_U2S_EL));


%sig2_ipp = (((CONST_F2^2/(CONST_F1^2-CONST_F2^2)).^2)*...
%             ( wrs2satdata(abv_mask,COL_U2S_SIG2L1MP)+...
%             wrs2satdata(abv_mask,COL_U2S_SIG2L2MP) ) + ...
%             GIVE_SIG2_L1L2REC + GIVE_SIG2_L1L2SAT)...
%             ./ obliquity2(wrs2satdata(abv_mask,COL_U2S_EL));
            
%sig2_ipp=diag(sig2_ipp);   
            
            
            
            
sig2_ipp = diag((((CONST_F2^2/(CONST_F1^2-CONST_F2^2)).^2)*...
             ( wrs2satdata(abv_mask,COL_U2S_SIG2L1MP)+...
             wrs2satdata(abv_mask,COL_U2S_SIG2L2MP) )) ...
             .*inv_obl2 );
%fill in the off-diagonal elements
obl_mat=sqrt(inv_obl2*inv_obl2');
sat_idx=wrs2satdata(abv_mask,COL_U2S_PRN);
wrs_s=unique(wrs_idx);
for idx =1:length(wrs_s);
     elem=find(wrs_idx==wrs_s(idx));
     sig2_ipp(elem,elem)=sig2_ipp(elem,elem) + ...
                  GIVE_SIG2_L1L2REC*obl_mat(elem,elem);
end
prn_s=unique(sat_idx);
for idx =1:length(prn_s);
     elem=find(sat_idx==prn_s(idx));
     sig2_ipp(elem,elem)=sig2_ipp(elem,elem) + ...
                  GIVE_SIG2_L1L2SAT*obl_mat(elem,elem);
end
% sig2_ipp = ones(length(abv_mask), 1) .* .0001;

if TRUTH_FLAG
  Ivpp = wrs2satdata(abv_mask,COL_U2S_IVPP);
end
%initialize return values
igpdata(:,COL_IGP_MINMON)=zeros(n_igp,1);
givei=repmat(MOPS_GIVEI_NM,n_igp,1);
upm_givei=givei-2;
betai=repmat(MOPS_GIVEI_NM,n_igp,1);

%loop over each Ionospheric Grid Point
for igp=1:n_igp
  %calculate the covariance matrix, radius, relative centroid metric and get
  %the indices for the IPPs used in the planar fit
  beta=0;
  if TRUTH_FLAG
      
      [cov, r, rcm, idx, delay, chi2] = igp_plane2(xyz_igp(igp,:), ...
                                            igp_en_hat(igp,:), xyz_ipp, ...
                                            sig2_ipp, Ivpp);
      if delay < 0
        delay = 0;
      end
      igpdata(igp,COL_IGP_DELAY) = delay;
  else
      [cov, r, rcm, idx,beta] = igp_plane(xyz_igp(igp,:), igp_en_hat(igp,:),...
                                     xyz_ipp, sig2_ipp);
  end
                                     
  % set rcm to the ipp_spread metric instead of relative centroid if
  % desired.  
  if (IPP_SPREAD_FLAG) 
      a_spd = 7;  % a_spd, r_spd tuning.
      r_spd = 424 * 1000;  % not used at the moment.
      diffs = xyz_ipp(idx,:) - ones(size(idx, 1),1) * xyz_igp(igp,:);
      dx = diffs * igpdata(igp, COL_IGP_EHAT)';
      dy = diffs * igpdata(igp, COL_IGP_NHAT)';
      
      rcm = obtain_ipp_spread(a_spd, r_spd, dx, dy, sig2_ipp(idx), GIVE_SIG2_DECORR);
  end;
                                 
  n_ipp=length(idx);
  %make sure we have the minimum number of points

  if(n_ipp>=GIVE_N_MIN)
	%check for at least 2 WRSs with 2 IPPs
    if(checkfor2(wrs_idx(idx)))
	  % flag that this meets minimum monitoring requirements
	  igpdata(igp,COL_IGP_MINMON)=1;
      
      if RTR_FLAG 
            %Rirreg2=RTR_MULT(n_ipp-3)*chi2;
            %Compute Raytheon Rirreg
            
            Rirreg2=chi2*GIVE_R2_DYN_IRREG(n_ipp-3);
            
      else
            Rirreg2=GIVE_R2_IRREG(n_ipp-3);
      end
      
      Rirreg2=max(1,Rirreg2);
      
      %calculate the maximum covariance projection at the 4 corners
      sig2_maxcorner=Rirreg2*max([...
             igp_corner_den(1,:,igp)*cov*igp_corner_den(1,:,igp)',...
             igp_corner_den(2,:,igp)*cov*igp_corner_den(2,:,igp)',...
             igp_corner_den(3,:,igp)*cov*igp_corner_den(3,:,igp)',...
             igp_corner_den(4,:,igp)*cov*igp_corner_den(4,:,igp)']);

      %calculate the max of the worst case decorr and undersampled threat
      
      
      sig2_threat=max([Rirreg2*GIVE_SIG2_DECORR ...
                    undersamp_threat(r/1000, rcm, GIVE_UNDERSAMP_PTS)]);
      %sig2_threat=Rirreg2*GIVE_SIG2_DECORR;
      
      %assemble GIVE and apply the ceiling function
      sig2_give=min([GIVE_KHMI_KA*(sig2_maxcorner + sig2_threat +...
                       GIVE_SIG2_ROT0), ...
                       sig2_max_vtec(mod(t + igp_lon(igp)*...
                       CONST_SEC_PER_DAY/360,CONST_SEC_PER_DAY), 0, ...
                       igp_mag_lat(igp), GIVE_KIGD_FAC, GIVE_CEILING_P, ...
                       GIVE_CEILING_T, GIVE_CEILING_LAT)]);

      %discretize the GIVE value
      if(sig2_give <= MOPS_SIG2_GIVE(MOPS_GIVEI_NM-1))
        temp=find(floor(MOPS_SIG2_GIVE/sig2_give));
        givei(igp)=temp(1);
        betai(igp)=beta;
        % UPM GIVE
        temp=find(floor(MOPS_SIG2_GIVE/(sig2_maxcorner + GIVE_SIG2_DECORR)));
        upm_givei(igp)=temp(1);
        
        %storm detector
        if TRUTH_FLAG
            if chi2 >= GIVE_CHI2_THRESH(n_ipp-3)
              TRIP_COUNT = TRIP_COUNT + 1;
              givei(igp) = MOPS_GIVEI_NM-1;
          end
        end
      end   
    end
  end
end

below_floor=find(givei < GIVE_GIVEI_FLOOR);
if(~isempty(below_floor))
  givei(below_floor)=GIVE_GIVEI_FLOOR;
end

igpdata(:,COL_IGP_GIVEI) = givei;
igpdata(:,COL_IGP_UPMGIVEI) = upm_givei;
igpdata(:,COL_IGP_BETA)=betai;





