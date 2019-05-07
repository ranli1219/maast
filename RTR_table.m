function table=RTR_table(N,phmi,K)
%*************************************************************************
%*     Copyright c 2007 The board of trustees of the Leland Stanford     *
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
%This function computes chi-square lower bound.  
%It solves the equation for all n up to N:
%f(chi2_lower_bound,n,K)=phmi
%K is usually 5.592
%phmi is the PHMI we want to protect against
%n is the number of degrees of freedom

%Juan Blanch, Sep 2007

alpha=zeros(N,1);
for k=1:N   
    alpha(k) = (finv(phmi,k,1)*k)*(K^2);  
end

table=alpha;
