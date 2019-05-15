function satdata = af_geoadd2(satdata, ~, wrs2satdata, ~, dual_freq)

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

% Created 2003 Mar19 By Todd Walter
% Added Antenna Bias 2004 May21 by Todd Walter
% Updated to Release 8/9 2009 Nov 20 by Todd Walter
% Updated to Release 3A1 2013 Aug 30 by Todd Walter

global GEO2_SIG2_ROT100 GEO2_SIG2_WREL1L2 
global GEO2_MIN_NWRS GEO2_SIG2_WRECLK 
global GEO2_UDREI_FLOOR GEO2_CEIL_EFF_GIVE
global GEO2_K1 GEO2_K2 GEO2_K3 GEO2_K4
global GEO2_DW_ELEV GEO2_DW1 GEO2_DW2 GEO2_DW3  GEO_CNMP_UNOB_BIAS_INIT
global GEO_CNMP_UNOB_BIAS GEO_CNMP_UNOB_BIAS_INIT_CNMP_VAL
global MOPS_UDREI_NM  
global COL_SAT_UDREI COL_SAT_COV   ...
         COL_SAT_MINMON COL_U2S_GXYZB COL_U2S_EL COL_U2S_IVPP ...
        COL_U2S_SIG2TRP COL_U2S_SIG2L1MP COL_U2S_SIG2L2MP 
global CONST_F1 CONST_F5

nsat = size(satdata,1);
% nwrs = size(wrsdata,1);
nlos = size(wrs2satdata,1);

sig2_trop = wrs2satdata(:,COL_U2S_SIG2TRP);

if(dual_freq)

    % assign L5 values to be equal to L1 values
    sig2_mon = sig2_trop + wrs2satdata(:,COL_U2S_SIG2L1MP)*...
               (CONST_F1^4 + CONST_F5^4)/((CONST_F1^2 - CONST_F5^2)^2) +...
                       GEO2_SIG2_WRECLK + GEO2_SIG2_WREL1L2; 
    %temporary  assign L2 values to be equal to L1 values
    %beacause COL_U2S_SIG2L2MP location is used to  store sigma_UIVE_UPM
    wrs2satdata(:,COL_U2S_SIG2L2MP) = wrs2satdata(:,COL_U2S_SIG2L1MP);  
    F2 = zeros(size((wrs2satdata(:,COL_U2S_EL))));
    
else
    F2 = obliquity2(wrs2satdata(:,COL_U2S_EL));
    %interpolated UPM_UIVE stored in COL_U2S_SIG2L2MP
    sig2_mon = sig2_trop + wrs2satdata(:,COL_U2S_SIG2L1MP) +...
                F2.*(GEO2_SIG2_ROT100 + wrs2satdata(:,COL_U2S_SIG2L2MP)) + ...
                GEO2_SIG2_WRECLK + GEO2_SIG2_WREL1L2;
end
mu_ant = GEO2_K1*(GEO2_K2 + GEO2_K3*exp(-wrs2satdata(:,COL_U2S_EL)*180/(pi*GEO2_K4)));

geo_unob_bias = GEO_CNMP_UNOB_BIAS_INIT*ones(nlos,1);
idx = find(wrs2satdata(:,COL_U2S_SIG2L1MP) <= ...
                     GEO_CNMP_UNOB_BIAS_INIT_CNMP_VAL);
if ~isempty(idx)
    geo_unob_bias(idx) = GEO_CNMP_UNOB_BIAS;
end
        
satdata(:,COL_SAT_MINMON)=zeros(nsat,1);
sat_udrei = repmat(MOPS_UDREI_NM,nsat,1);
for isat = 1:nsat
    % extract los indices for each satellite
    wrs4satidx = isat:nsat:nlos;

    % find satellites satisfying minimum visibility
    % not vis if el is NaN.
    idxvis = find(~isnan(wrs2satdata(wrs4satidx,COL_U2S_EL)) & ...
                   (wrs2satdata(wrs4satidx,COL_U2S_SIG2L2MP) > 0)); 
    
    % calculate udre for those that satisfy min visibility
    if (length(idxvis)>= GEO2_MIN_NWRS)
        idxvis = wrs4satidx(idxvis);    % refer back to wrs2satdata indexing

		satdata(isat,COL_SAT_MINMON)=1;
        sig2_udre_cp = (4.5/3.29).^2;
        
        idx_good_give = find(wrs2satdata(idxvis,COL_U2S_IVPP) > 0);

        min_uire2_sp = min([min(F2(idxvis(idx_good_give)).* ...
                      wrs2satdata(idxvis(idx_good_give),COL_U2S_IVPP)) ...
                            GEO2_CEIL_EFF_GIVE]);
        
        % Calculate covariance matrix and psi_max
        G = -wrs2satdata(idxvis,COL_U2S_GXYZB);
        Pmt28 = reshape(satdata(isat,COL_SAT_COV)',4,4);
        psi_max = mt28geo_cov2(G,sig2_mon(idxvis),chol(Pmt28));
        
        % Calculate deweighting factor
        alpha_dw=ones(size(idxvis))';
        lowev=find(wrs2satdata(idxvis, COL_U2S_EL) <= GEO2_DW_ELEV);
        if(~isempty(lowev))
          alpha_dw(lowev)= sqrt(1./(GEO2_DW1*wrs2satdata(idxvis(lowev), ...
                    COL_U2S_EL)*180/pi + GEO2_DW2) + GEO2_DW3);  
        end                
        % Calculate UDRE by phmi equation
        sat_udrei(isat) = geo_phmi_eq2(sig2_udre_cp, sig2_mon(idxvis), ...
                sig2_trop(idxvis), mu_ant(idxvis), alpha_dw, psi_max, ...
                min_uire2_sp, geo_unob_bias(idxvis));

        if(sat_udrei(isat) < GEO2_UDREI_FLOOR)
          sat_udrei(isat) = GEO2_UDREI_FLOOR;
        end
    end

end   
satdata(:,COL_SAT_UDREI) = sat_udrei;









