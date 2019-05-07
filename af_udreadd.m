function satdata = af_udreadd(satdata,wrsdata,wrs2satdata,do_mt28)

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
% fixed square root in low elevation deweighting 2009 Nov 23 Todd Walter

global UDRE_MIN_NWRS UDRE_MIN_NWRSCONUS UDRE_SIG2_WRECLK 
global UDREI_FLOOR UDREI_FLOOR_NOMINCONUS
global UDRE_DW_ELEV UDRE_DW1 UDRE_DW2 UDRE_DW3
global UDRE_K1 UDRE_K2 UDRE_K3 UDRE_K4 UDRE_K5 UDRE_K6 UDRE_K7 UDRE_K8
global CONST_H_IONO CONST_F1 CONST_F2
global MOPS_UDREI_NM  
global COL_USR_INBND COL_SAT_UDREI COL_SAT_COV COL_SAT_XYZ COL_SAT_XYZDOT ...
        COL_SAT_SCALEF COL_SAT_MINMON COL_U2S_GXYZB COL_U2S_EL ...
        COL_U2S_SIG2TRP COL_U2S_SIG2L1MP COL_U2S_SIG2L2MP COL_SAT_SIG2CP

nsat = size(satdata,1);
nwrs = size(wrsdata,1);
nlos = size(wrs2satdata,1);

sig2_wreclk = repmat(UDRE_SIG2_WRECLK,nlos,1);

a1 = CONST_F1^2/(CONST_F1^2-CONST_F2^2);
a2 = CONST_F2^2/(CONST_F1^2-CONST_F2^2);
sig2_mon = wrs2satdata(:,COL_U2S_SIG2TRP) + ...
            a1^2.*wrs2satdata(:,COL_U2S_SIG2L1MP) + ...
            a2^2.*wrs2satdata(:,COL_U2S_SIG2L2MP) + sig2_wreclk;
   
el = wrs2satdata(:,COL_U2S_EL)*180/pi;
L1_ave = UDRE_K1*(UDRE_K2 + UDRE_K3*exp(-el/UDRE_K4));
L2_ave = UDRE_K5*(UDRE_K6 + UDRE_K7*exp(-el/UDRE_K8));
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
    
    % find number in conus
    idxconus = find(wrsdata(:,COL_USR_INBND)); 
    nconus=length(findcommon(idxvis,idxconus));

    % calculate udre for those that satisfy min visibility
    if (length(idxvis)>= UDRE_MIN_NWRS)
        idxvis = wrs4satidx(idxvis);    % refer back to wrs2satdata indexing

		satdata(isat,COL_SAT_MINMON)=1;
        sig2_udre_cp = (-0.0134*length(idxvis) + 0.9521).^2;

        % Calculate MT28 normalized covariance matrix
        if do_mt28,
            G = -wrs2satdata(idxvis,COL_U2S_GXYZB);
            [dCov,scalef,maxrat] = mt28cov(satdata(isat,COL_SAT_XYZ),...
                            satdata(isat,COL_SAT_XYZDOT),G,sig2_mon(idxvis));
            satdata(isat,COL_SAT_COV)=dCov(:)';
            satdata(isat,COL_SAT_SCALEF)=scalef;
        else 
            maxrat = 1;
        end

        % Calculate deweighting factor
        alpha=ones(size(idxvis))';
        lowev=find(wrs2satdata(idxvis, COL_U2S_EL) <= UDRE_DW_ELEV);
        if(~isempty(lowev))
          alpha(lowev)= sqrt(1./(UDRE_DW1*wrs2satdata(idxvis(lowev), ...
                          COL_U2S_EL)*180/pi + UDRE_DW2) + UDRE_DW3);
        end
        % Calculate UDRE by phmi equation
        sat_udrei(isat) = phmi_eq(sig2_udre_cp, sig2_mon(idxvis), ...
                                  mu_ant(idxvis), alpha, maxrat);
        if(nconus >= UDRE_MIN_NWRSCONUS)
          udrei_floor = UDREI_FLOOR;
        else
          udrei_floor = UDREI_FLOOR_NOMINCONUS;
        end
        if(sat_udrei(isat) < udrei_floor)
          sat_udrei(isat) = udrei_floor;
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






