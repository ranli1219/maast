function satdata = af_udregeoadd(satdata,wrsdata,wrs2satdata, do_mt28, dual_freq)

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

% Created 2003 Mar19 By Todd Walter
% Added Antenna Bias 2004 May21 by Todd Walter

global GEO_SIG2_CNMPFLOOR  GEO_SIG2_ROT100 GEO_SIG2_WREL1L2 GEO_SIG2_WRECLK
global GEO_MIN_NWRS GEO_MIN_NWRSCONUS GEO_SIG2_WRECLK 
global GEO_UDREI_FLOOR GEO_UDREI_FLOOR_NOMINCONUS
global GEO_K1 GEO_K2 GEO_K3 GEO_K4
global MOPS_UDREI_NM  
global COL_USR_INBND COL_SAT_UDREI COL_SAT_COV COL_SAT_XYZ COL_SAT_XYZDOT ...
        COL_SAT_SCALEF COL_SAT_MINMON COL_U2S_GXYZB COL_U2S_EL ...
        COL_U2S_SIG2TRP COL_U2S_SIG2L1MP COL_U2S_SIG2L2MP COL_SAT_SIG2CP
global CONST_F1 CONST_F5

nsat = size(satdata,1);
nwrs = size(wrsdata,1);
nlos = size(wrs2satdata,1);

if(dual_freq)

    %temporary  assign L2 values to be equal to L1 values
    wrs2satdata(:,COL_U2S_SIG2L2MP) = GEO_SIG2_CNMPFLOOR;
    sig2_mon = wrs2satdata(:,COL_U2S_SIG2TRP) + GEO_SIG2_CNMPFLOOR*...
               (CONST_F1^4 + CONST_F5^4)/((CONST_F1^2 - CONST_F5^2)^2) +...
                       GEO_SIG2_WRECLK + GEO_SIG2_WREL1L2; 
else
    F2 = obliquity2(wrs2satdata(:,COL_U2S_EL));
    %interpolated UPM_UIVE stored in COL_U2S_SIG2L2MP
    sig2_mon = wrs2satdata(:,COL_U2S_SIG2TRP) + GEO_SIG2_CNMPFLOOR +...
                F2.*(GEO_SIG2_ROT100 + wrs2satdata(:,COL_U2S_SIG2L2MP)) + ...
                GEO_SIG2_WRECLK + GEO_SIG2_WREL1L2;
end
mu_ant = GEO_K1*(GEO_K2 + GEO_K3*exp(-wrs2satdata(:,COL_U2S_EL)*180/(pi*GEO_K4)));

       
        
satdata(:,COL_SAT_MINMON)=zeros(nsat,1);
sat_udrei = repmat(MOPS_UDREI_NM,nsat,1);
for isat = 1:nsat,
    % extract los indices for each satellite
    wrs4satidx = isat:nsat:nlos;

    % find satellites not satisfying minimum visibility
    % not vis if el is NaN.
    idxvis = find(~isnan(wrs2satdata(wrs4satidx,COL_U2S_EL)) & ...
                   (wrs2satdata(wrs4satidx,COL_U2S_SIG2L2MP) > 0)); 
    
    % find number in conus
    idxconus = find(wrsdata(:,COL_USR_INBND)); 
    nconus=length(findcommon(idxvis,idxconus));

    % calculate udre for those that satisfy min visibility
    if (length(idxvis)>= GEO_MIN_NWRS)
        idxvis = wrs4satidx(idxvis);    % refer back to wrs2satdata indexing

		satdata(isat,COL_SAT_MINMON)=1;
%        sig2_udre_cp = (4.2/3.29).^2;
        sig2_udre_cp = (11/3.29).^2;
        
        % Calculate covariance matrix
        G = -wrs2satdata(idxvis,COL_U2S_GXYZB);
        Cov = inv(G'*diag(ones(size(sig2_mon(idxvis)))./sig2_mon(idxvis))*G );
        lambda2_wrs=diag(G*Cov*G');
        lambda2_wrs = lambda2_wrs/(max(lambda2_wrs)*1.1);
        
        % Calculate UDRE by phmi equation
        sat_udrei(isat) = geo_phmi_eq(sig2_udre_cp, sig2_mon(idxvis), ...
                           mu_ant(idxvis), lambda2_wrs);

        if(nconus >= GEO_MIN_NWRSCONUS)
          udrei_floor = GEO_UDREI_FLOOR;
        else
          udrei_floor = GEO_UDREI_FLOOR_NOMINCONUS;
        end
        if(sat_udrei(isat) < udrei_floor)
          sat_udrei(isat) = udrei_floor;
        end
    end

end   
satdata(:,COL_SAT_UDREI) = sat_udrei;

%if using MT 28 put in the identity matrix for XYZ and 0 for clock
if do_mt28
  a=eye(4);
  a(4,4)=0;
  a=a(:)';
  satdata(:,COL_SAT_COV)=repmat(a,nsat,1);
  satdata(:,COL_SAT_SCALEF)=repmat(0,nsat,1);
end   






