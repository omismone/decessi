clc; clearvars; close all;
raw_deaths = readtable("res\deceduti.csv");
raw_positives = readtable("res\positivi.csv"); 
% positives data dates start 22 days before deaths (see csv) so:
raw_positives(1:22,:) = [];


interval = 225:375;


deaths = table2array(raw_deaths(interval, 3)); % working only on weekly averages
positives = table2array(raw_positives(interval, 3));


% Plot
dates = table2array(raw_positives(interval,1));

% figure(1)
% plot(dates, positives);
% hold on
% plot(dates, deaths);
% legend("Positives","Deaths");



%% ROTTEN MODEL

% Determining the shift 
max_shift = 100;
possible_values = 1:max_shift;
cc_array = zeros(length(possible_values),1);
for i = possible_values
    i_shifted_positives = table2array(raw_positives(interval-i,3));
    cc = corrcoef(i_shifted_positives, deaths);
    cc_array(i) = cc(1,2);
end
shift = find(cc_array == max(cc_array));
shifted_positives = table2array(raw_positives(interval-shift, 3));


% Determining the gain
gain = lscov(shifted_positives,deaths);
deamplified_shifted_positives = shifted_positives.*gain;

ssr = (deaths - deamplified_shifted_positives)' * (deaths - deamplified_shifted_positives);

% Plot
figure(2)
plot(dates, deaths)
title("rotten model")
hold on
plot(dates, deamplified_shifted_positives)
legend("deaths","deamplified shifted positives")


% Checking all the possibilities
shift_range = linspace(1,224,224);
gain_range = linspace(0,0.1,224);
ssr_matrix = [zeros(size(shift_range))];

count_i = 1;
count_j = 1;
for i = shift_range
    for j = gain_range
        ssr_matrix(count_i,count_j) = calc_ssr(raw_positives, deaths, interval, i,j); %shift on rows
        count_j = count_j +1;
    end
    count_j = 1;
    count_i = count_i +1;
end
min = min(min(ssr_matrix));
[mountain_shift_index, mountain_gain_index] = find(ssr_matrix == min);
mountain_shift = shift_range(mountain_shift_index);
mountain_gain = gain_range(mountain_gain_index);
% Plot
figure(3)
[shift_grid, gain_grid] = meshgrid(shift_range,gain_range);
mesh(shift_grid, gain_grid, ssr_matrix') 
xlabel("shift")
ylabel("gain")
zlabel("ssr")
ylim([0,0.07])
zlim([-1,100000000])
xlim([0,80])
grid on
title("params mountain view")
% fprintf(sprintf("shift:" + shift + ", 'mountain' shift:"+mountain_shift+"\ngain:" + gain + ", 'mountain' gain:"+mountain_gain+"\n"));



% plotting separately
figure(4)
index = 58;
subplot(2,1,1)
plot(ssr_matrix(:,index))
xlim([0,224])
stringa = sprintf("shift (gain="+gain_range(index)+")");
xlabel(stringa)
ylabel("ssr")
subplot(2,1,2)
plot(gain_grid,ssr_matrix(index,:))
stringa = sprintf("gain (shift="+shift_range(index)+")");
xlabel(stringa)
ylabel("ssr")




%% 2 model






%% functions
function [ssr] = calc_ssr(raw_positives, deaths, interval, shift, gain)
    ssr = (deaths - (table2array(raw_positives(interval-shift,3)).*gain))' * (deaths - (table2array(raw_positives(interval-shift,3)).*gain));
end