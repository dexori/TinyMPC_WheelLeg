% SPDX-License-Identifier: MIT
function [L0,la,I,Y,fai_c] = balance_canshujisuan(theta1,theta2,AO,OO1,AB,BO1,l1,l2,ax,bx,Il1,Il2,ml1,ml2)
%%赋值
theta1=theta1*pi/180;
theta2=theta2*pi/180;
[L0,fai_c,A_GFJ]=balance_Leg(theta1,theta2,AO,OO1,AB,BO1,l1,l2);
%% 计算腿部质心
x1=ax*l1*cos(theta1);
x2=l1*cos(theta1)-bx*l2*cos(A_GFJ);

y1=ax*l1*sin(theta1);
y2=l1*sin(theta1)+bx*l2*sin(A_GFJ);

X=(x1*ml1+x2*ml2)/(ml1+ml2);
Y=(y1*ml1+y2*ml2)/(ml1+ml2);


%% 转动惯量--求各个杆
d1 = sqrt((X - x1)^2 + (Y - y1)^2) ;
d2 = sqrt((X - x2)^2 + (Y - y2)^2) ;

I = (Il1 + ml1 *d1^2) + (Il2 + ml2 * d2^2) ;%平行轴定理

%% 驱动轮到腿部质心距离
la=L0-(X-Y*tan(pi/2-fai_c))*sin(pi/2-fai_c)-Y/cos(pi/2-fai_c);

end
