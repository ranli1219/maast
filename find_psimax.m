function psi_max = find_psimax(P_normalized, L)

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
% See Appendix G of UDRE ADD A014-009L

% Created 2013 Aug 30 By Todd Walter

Sigma = [[1 0 0 0 ]; [0 1 0 0]; [0 0 1 0]; [0 0 0 -1]];
LI = inv(L);


A = LI'*P_normalized*LI;

B = LI'*Sigma*LI;

a = sort(eig(A));

b = sort(eig(B));

mu_min = (a(4) - a(1))/b(1);
mu_max = (a(4) - a(1))/b(4);
mu_mid = 0;

fmu_min = max(eig(A + mu_min*B));
fmu_max = max(eig(A + mu_max*B));
fmu_mid = a(4); 

cntr = 0;

while (abs(fmu_max - fmu_min) > 1e-4 && cntr < 15)
    mu_1 = (mu_min + mu_mid)/2;
    fmu_1 = max(eig(A + mu_1*B));
    
    mu_3 = (mu_mid + mu_max)/2;
    fmu_3 = max(eig(A + mu_3*B));
    
    if fmu_1 < fmu_mid && fmu_1 < fmu_3
        mu_mid = mu_1;
        fmu_mid = fmu_1;
        mu_max = mu_3;
        fmu_max = fmu_3;
    elseif  fmu_mid < fmu_1 && fmu_mid < fmu_3
        mu_min = mu_1;
        fmu_min = fmu_1;
        mu_max = mu_3;
        fmu_max = fmu_3;        
    else
        mu_min = mu_mid;
        fmu_min = fmu_mid;
        mu_mid = mu_3;
        fmu_mid = fmu_3;        
    end
    psi_max = min([fmu_min fmu_mid fmu_max]);
    cntr = cntr +1;
end
        










