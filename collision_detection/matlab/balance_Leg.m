% SPDX-License-Identifier: MIT
%% 串联腿连杆解算
function [L0,fai_c,A_GFJ,A_ABO,A_BAO,A_DEG] = balance_Leg(theta1,theta2,AO,OO1,AB,BO1,OF,GF)
    AO1=sqrt(AO.^2+OO1.^2-2*AO*OO1*cos(theta2));
    A_BAO1=acos((AB.^2+AO1.^2-BO1.^2)/(2*AB*AO1));
    A_ABO=acos((AB.^2+BO1.^2-AO1.^2)/(2*AB*BO1));
    A_O1AO=acos((AO1.^2+AO.^2-OO1.^2)/(2*AO1*AO));
    A_BAO=A_BAO1+A_O1AO;
    A_DEG=2*pi-theta2-A_BAO1-A_O1AO-A_ABO;
    A_GFJ=A_DEG-theta1;
    xc=OF*cos(theta1)-GF*cos(A_GFJ);
    yc=OF*sin(theta1)+GF*sin(A_GFJ);
    L0=sqrt(xc.^2+yc.^2);
    fai_c=atan2(yc,xc);
end
