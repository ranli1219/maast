function satdata = af_udreadd2(satdata,wrsdata,wrs2satdata,do_mt28)

%*************************************************************************
%*     Copyright c 2009 The board of trustees of the Leland Stanford     *
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
% Antenna bias added By Todd Walter May 21, 2004
% Updated to Release 8/9 2009 Nov 23 Todd Walter

global UDRE2_MIN_NWRS UDRE2_SIG2_WRECLK 
global UDRE2_UDREI_FLOOR
global UDRE2_DW_ELEV UDRE2_DW1 UDRE2_DW2 UDRE2_DW3
global UDRE2_K1 UDRE2_K2 UDRE2_K3 UDRE2_K4 UDRE2_K5 UDRE2_K6 UDRE2_K7 UDRE2_K8
global CONST_F1 CONST_F2
global MOPS_UDREI_NM  
global COL_SAT_UDREI COL_SAT_COV COL_SAT_XYZ COL_SAT_XYZDOT ...
        COL_SAT_SCALEF COL_SAT_MINMON COL_U2S_GXYZB COL_U2S_EL ...
        COL_U2S_SIG2TRP COL_U2S_SIG2L1MP COL_U2S_SIG2L2MP 

nsat = size(satdata,1);
nwrs = size(wrsdata,1);
nlos = size(wrs2satdata,1);

sig2_wreclk = repmat(UDRE2_SIG2_WRECLK,nlos,1);

a1 = CONST_F1^2/(CONST_F1^2-CONST_F2^2);
a2 = CONST_F2^2/(CONST_F1^2-CONST_F2^2);
sig2_trop = wrs2satdata(:,COL_U2S_SIG2TRP);
sig2_mon =  sig2_trop + ...
            a1^2.*wrs2satdata(:,COL_U2S_SIG2L1MP) + ...
            a2^2.*wrs2satdata(:,COL_U2S_SIG2L2MP) + sig2_wreclk;
   
el = wrs2satdata(:,COL_U2S_EL)*180/pi;
L1_ave = UDRE2_K1*(UDRE2_K2 + UDRE2_K3*exp(-el/UDRE2_K4));
L2_ave = UDRE2_K5*(UDRE2_K6 + UDRE2_K7*exp(-el/UDRE2_K8));
mu_ant = a1*L1_ave + a2*L2_ave;

satdata(:,COL_SAT_MINMON)=zeros(nsat,1);
sat_udrei = repmat(MOPS_UDREI_NM,nsat,1);
cnt=0;
for isat = 1:nsat,
    % extract los indices for each satellite
    wrs4satidx = isat:nsat:nlos;

    % find satellites not satisfying minimum visibility
    % not vis if el is NaN.
    idxvis = find(~isnan(wrs2satdata(wrs4satidx,COL_U2S_EL))); 
    

    % calculate udre for those that satisfy min visibility
    if (length(idxvis)>= UDRE2_MIN_NWRS)
        idxvis = wrs4satidx(idxvis);    % refer back to wrs2satdata indexing

		satdata(isat,COL_SAT_MINMON)=1;
        sig2_udre_cp = (-0.0134*length(idxvis) + 0.9521).^2;

        % Calculate MT28 normalized covariance matrix
        G = -wrs2satdata(idxvis,COL_U2S_GXYZB);
        [dCov,scalef,maxrat] = mt28cov2(satdata(isat,COL_SAT_XYZ),...
                          satdata(isat,COL_SAT_XYZDOT),G,sig2_mon(idxvis));
        satdata(isat,COL_SAT_COV)=dCov(:)';
        satdata(isat,COL_SAT_SCALEF)=scalef;
        
        % Calculate deweighting factor
        alpha_dw=ones(size(idxvis))';
        lowev=find(wrs2satdata(idxvis, COL_U2S_EL) <= UDRE2_DW_ELEV);
        if(~isempty(lowev))
          alpha_dw(lowev)= sqrt(1./(UDRE2_DW1*wrs2satdata(idxvis(lowev),...
                      COL_U2S_EL)*180/pi + UDRE2_DW2) + UDRE2_DW3);
        end
        % Calculate UDRE by phmi equation
        sat_udrei(isat) = phmi_eq2(sig2_udre_cp, sig2_mon(idxvis), ...
                      sig2_trop(idxvis), mu_ant(idxvis), alpha_dw, maxrat);
        
        if(sat_udrei(isat) < UDRE2_UDREI_FLOOR)
          sat_udrei(isat) = UDRE2_UDREI_FLOOR;
        end
    end

    if sat_udrei(isat)<MOPS_UDREI_NM,
        fprintf('o');
    else
        fprintf('.');
    end
    if mod(cnt,72)==71
        fprintf('\n');
    end
    cnt=cnt+1;
end   
satdata(:,COL_SAT_UDREI) = sat_udrei;
fprintf('\n');






