%POWER LAW
clear all
station = readtable('20070714093540r.csv');

station = removevars(station, {'Sample_'});
station = removevars(station, {'Date_Time'});
station = removevars(station, {'Frequency_MHz_'});
station = removevars(station, {'ProfileType'});
station = table2array(station);
station(station==0) = NaN;

[R, C] = size(station);

average_values = nanmean(station);
average_values = average_values';

cell_location = average_values(4:7:end);
%cell_location = cell_location .* (-1);

Ve = average_values(5:7:end);
Vn = average_values(6:7:end);
Vu = average_values(7:7:end);
Vd = average_values(8:7:end);
Vmag = sqrt(Ve.^2 + Vn.^2 + Vu.^2);

y2 = cell_location;
x2 = Vmag;

last_point = cell_location(end);

% figure(1);
% plot(flipud(x2), y2, 'ro')
% y2 = y2.* -1;
% y2 = sort(y2, 'descend');
% y2 = flipud(y2);
% x2 = flipud(x2);

% y2(1) = [];
% x2(1) = [];
% 
% y2(2) = [];
% x2(2) = [];

Q = 0;
x2 = x2'
x2 = horzcat(Q, x2);
x2(end+1) = 0;
x2 = x2';

y2 = flipud(y2);
y2(end+1) = 0;
y2 = flipud(y2);
y2(end+1) = 3.13;
y2 = y2 ./ max(y2);

% y2(54) = [];
% x2(54) = [];

z_bed = y2(end-7:end-1);
u_bed = x2(end-7:end-1);
z_ice = y2(2:7);
u_ice = x2(2:7);

z_bed_zero = y2(end-7:end);
u_bed_zero = x2(end-7:end);
z_ice_zero = y2(1:7);
u_ice_zero = x2(1:7);

z_bed = flipud(z_bed);
u_bed = flipud(u_bed);
z_ice = flipud(z_ice);
u_ice = flipud(u_ice);

z_bed_zero = flipud(z_bed_zero);
u_bed_zero = flipud(u_bed_zero);
z_ice_zero = flipud(z_ice_zero);
u_ice_zero = flipud(u_ice_zero);

c_bed = polyfit(z_bed,u_bed,1);
disp(['Equation for bed layer y = ' num2str(c_bed(1)) '*x + ' num2str(c_bed(2))])
y_est_bed = polyval(c_bed,z_bed);
slope_bed = (y_est_bed(2) - y_est_bed(1)) / (z_bed(2) - z_bed(1));

c_ice = polyfit(z_ice,u_ice,1);
disp(['Equation for ice layer y = ' num2str(c_ice(1)) '*x + ' num2str(c_ice(2))])
y_est_ice = polyval(c_ice,z_ice);
slope_ice = (y_est_ice(2) - y_est_ice(1)) / (z_ice(2) - z_ice(1));

c_bed_zero = polyfit(z_bed_zero,u_bed_zero,1);
disp(['Equation for bed layer with zero point y = ' num2str(c_bed_zero(1)) '*x + ' num2str(c_bed_zero(2))])
y_est_bed_zero = polyval(c_bed_zero,z_bed_zero);
% slope_bed_zero = (y_est_bed_zero(2) - y_est_bed_zero(1)) / (z_bed_zero(2) - z_bed_zero(1));

c_ice_zero = polyfit(z_ice_zero,u_ice_zero,1);
disp(['Equation for ice layer with zero point y = ' num2str(c_ice_zero(1)) '*x + ' num2str(c_ice_zero(2))])
y_est_ice_zero = polyval(c_ice_zero,z_ice_zero);
% slope_ice_zero = (y_est_ice_zero(2) - y_est_ice_zero(1)) / (z_ice_zero(2) - z_ice_zero(1));




%POWER LAW
k_s = 0.55;
slope_bed_zero = -1.3;
slope_ice_zero = 1.3;

ice_term = (y2 ./ max(y2)) .^ (1/slope_ice_zero);
bed_term = (1 - y2 ./ max(y2)) .^ (1/-slope_bed_zero);

u_power = k_s .* ice_term .* bed_term;


Vmag(end+1) = 0;
Vmag = flipud(Vmag);
Vmag(end+1) = 0;
Vmag = flipud(Vmag);

figure(1);
plot(flipud(x2), y2, 'ro')
hold on
plot(y_est_bed, 1-z_bed, 'k--','LineWidth',2)
hold on
plot(y_est_ice, 1-z_ice, 'g--','LineWidth',2)
hold on
plot(y_est_bed_zero, 1-z_bed_zero, 'k--','LineWidth',2)
hold on
plot(y_est_ice_zero, 1-z_ice_zero, 'g--','LineWidth',2)
hold on
plot(flipud(u_power), y2 ./ max(y2), 'b-', 'LineWidth',2)
legend('original data', 'bed slope', 'iceslope', 'bed slope with 0', 'ice slope with 0', 'power law', 'Location','west');

% figure(2);
% plot(u_power, y2 ./ max(y2), 'b-', 'LineWidth',2)

% figure(2);
% plot(u_power, y2 ./ max(y2), 'b-')
% 
% %Extra point at Z = 0.6
% extra_point_u1 = y_est(1) - (slope_bed * (fit_Z(1) - 2.5));

%figure(2);
% plot(z_bed, u_bed, 'ko')
% xlabel('z (m)')
% ylabel('Velocity (m/s)')
% hold on
% plot(fit_Z,y_est,'r--','LineWidth',2)
% xlim([0 4])
% ylim([0.05 0.5])
% set(gca,'XTick',(0:2:4))
% xline(first_20_percent, 'r--', 'LineWidth', 1);
% text(first_20_percent,0.04,'0.2 z/H','FontSize',10)
% hold on
% xline(first_50_percent, 'r--', 'LineWidth', 1);
% text(first_50_percent,0.04,'0.5 z/H','FontSize',10)
% ylim([0.1 0.25])
% hold off