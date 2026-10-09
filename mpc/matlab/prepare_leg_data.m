function leg = prepare_leg_data(p)
%PREPARE_LEG_DATA 读取外部 CAD 表，计算每个机构角的腿部等效参数。
%
% CSV 第一行必须是：
% dist_mm,com_y_mm,com_x_mm,Ixx_kgmm2

fid = fopen(p.mass_property_csv,'r');
if fid < 0
    error('Cannot open CSV: %s',p.mass_property_csv);
end
header = strtrim(fgetl(fid));
fclose(fid);
required_header = 'dist_mm,com_y_mm,com_x_mm,Ixx_kgmm2';
if ~strcmp(header,required_header)
    error('CSV header must be: %s',required_header);
end

table_data = dlmread(p.mass_property_csv,',',1,0);
if size(table_data,2) ~= 4 || any(~isfinite(table_data(:)))
    error('CSV must contain four numeric columns.');
end
if any(diff(table_data(:,1)) <= 0)
    error('CSV dist_mm must be strictly increasing.');
end

dist = table_data(:,1)*1e-3;
com_y = -table_data(:,2)*1e-3;  % 沿用原工程的 CAD 坐标正负号
com_x = table_data(:,3)*1e-3;
inertia = table_data(:,4)*1e-6;

angles = p.angle_range_deg(:);
for i = 1:numel(angles)
    angle1 = angles(i);
    angle4 = 180-2*angle1;
    L0 = fivebar_length(angle1,angle4,p);
    if L0 < dist(1) || L0 > dist(end)
        error('Leg length %.6f m is outside the CSV range.',L0);
    end

    x = interp1(dist,com_x,L0,'pchip');
    y = interp1(dist,com_y,L0,'pchip');
    leg(i).angle_deg = angle1; %#ok<AGROW>
    leg(i).length = L0;
    leg(i).inertia = interp1(dist,inertia,L0,'pchip');
    leg(i).body_to_com = hypot(x,y);
    leg(i).wheel_to_com = hypot(x,L0-y);
    leg(i).body_offset = atan2(x,y);
    leg(i).wheel_offset = atan2(x,L0-y);
end

lengths = [leg.length];
fprintf('  CSV range %.6f .. %.6f m\n',dist(1),dist(end));
fprintf('  Leg range %.6f .. %.6f m\n',min(lengths),max(lengths));
end

function L0 = fivebar_length(theta1_deg,theta2_deg,p)
% 与原工程 balance_Leg.m 相同的五连杆几何。
theta1 = theta1_deg*pi/180;
theta2 = theta2_deg*pi/180;
AO1 = sqrt(p.AO^2+p.OO1^2-2*p.AO*p.OO1*cos(theta2));
A_BAO1 = acos((p.AB^2+AO1^2-p.BO1^2)/(2*p.AB*AO1));
A_ABO = acos((p.AB^2+p.BO1^2-AO1^2)/(2*p.AB*p.BO1));
A_O1AO = acos((AO1^2+p.AO^2-p.OO1^2)/(2*AO1*p.AO));
A_DEG = 2*pi-theta2-A_BAO1-A_O1AO-A_ABO;
A_GFJ = A_DEG-theta1;
xc = p.l1*cos(theta1)-p.l2*cos(A_GFJ);
yc = p.l1*sin(theta1)+p.l2*sin(A_GFJ);
L0 = hypot(xc,yc);
end
