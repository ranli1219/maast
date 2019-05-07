function [p,a0,a1,a2,a3]=char_poly(l,A)
%*************************************************************************
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
n=length(l);
a=A(1,1);
b=A(1,2);
c=A(1,3);
d=A(1,4);
e=A(2,2);
f=A(2,3);
g=A(2,4);
h=A(3,3);
i=A(3,4);
j=A(4,4);
a0=-2*b*g*c*i+2*b*g*d*h+(b*i)^2-e*(d^2)*h+(c*g)^2-e*(c^2)*j+(d*f)^2+a*e*h*j+2*c*e*d*i...
   -2*c*g*d*f+2*b*f*c*j-a*(g^2)*h-2*b*f*d*i-a*e*(i^2)-(b^2)*h*j+2*a*f*g*i-a*(f^2)*j;
a1=-2*f*g*i-a*e*j+e*(c^2)-2*c*d*i+(d^2)*h+a*(f^2)-2*b*g*d+(b^2)*h+(b^2)*j+(g^2)*h+e*(i^2)...
   -a*h*j+a*(g^2)+a*(i^2)+e*(d^2)-e*h*j+(f^2)*j-a*e*h-2*b*f*c+(c^2)*j;
a2=a*h+h*j+a*j+e*h-b^2+e*j-c^2-f^2-d^2+a*e-i^2-g^2;
a3=-h-j-a-e;
p=(l^4)+a3*(l^3)+a2*(l^2)+a1*l+a0*eye(n);
