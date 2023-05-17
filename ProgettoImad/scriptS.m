clc; clearvars; close all;
raw_deaths = readtable('iss_bydate_italia_deceduti.csv');
raw_positives = readtable('iss_bydate_italia_positivi.csv'); 
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



%% The Rotten Model

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
legend("deaths","estimated deaths")


% Checking all the possibilities
shift_range = linspace(1,224,224);
gain_range = linspace(0,0.1,224);
ssr_matrix = [zeros(size(shift_range))];

count_i = 1;
count_j = 1;
for i = shift_range
    for j = gain_range
        estimation = (table2array(raw_positives(interval-i,3)).*j); %shift is on rows
        ssr_matrix(count_i,count_j) = getSsr(deaths, estimation); 
        count_j = count_j +1;
    end
    count_j = 1;
    count_i = count_i +1;
end
minimum = min(ssr_matrix(:));
[grid_shift_index, grid_gain_index] = find(ssr_matrix == minimum);
grid_shift = shift_range(grid_shift_index);
grid_gain = gain_range(grid_gain_index);
% Plot
figure(3)
[shift_grid, gain_grid] = meshgrid(shift_range,gain_range);
mesh(shift_grid, gain_grid, ssr_matrix') 
xlabel("shift")
ylabel("gain")
zlabel("ssr")
xlim([shift_range(1),shift_range(length(shift_range))])
ylim([gain_range(1),gain_range(length(gain_range))])
zlim([minimum,max(ssr_matrix(:))])
grid on
title("rotten model parameter grid")
% fprintf(sprintf("shift:" + shift + ", 'grid' shift:"+grid_shift+"\ngain:" + gain + ", 'grid' gain:"+grid_gain+"\n"));



% plotting separately
% figure(4)
% index = 58;
% subplot(2,1,1)
% plot(ssr_matrix(:,index))
% xlim([0,224])
% stringa = sprintf("shift (gain="+gain_range(index)+")");
% xlabel(stringa)
% ylabel("ssr")
% subplot(2,1,2)
% plot(gain_grid,ssr_matrix(index,:))
% stringa = sprintf("gain (shift="+shift_range(index)+")");
% xlabel(stringa)
% ylabel("ssr")



%% 2 The Good Model
 
% Il ritardo è fissato, cambiarne il valore a mano (forma una curva ad U con ssr minimo a D = 13).
D = 13;
lambda_range = linspace(1,3,100);
gain_range = linspace(0,0.1,100);
ssr_matrix = [zeros(size(lambda_range))];

count_i = 1;
count_j = 1;
u = table2array(raw_positives(:,3));
for i = lambda_range
    for j = gain_range
        estimation = getGoodModelYCap(u,interval,D,j,i);
        ssr_matrix(count_i,count_j) = getSsr(deaths, estimation); 
        count_j = count_j +1;
    end
    count_j = 1;
    count_i = count_i +1;
end
minimum = min(ssr_matrix(:));
[grid_lambda_index, grid_gain_index] = find(ssr_matrix == minimum);
grid_lambda = lambda_range(grid_lambda_index);
grid_gain = gain_range(grid_gain_index);
% Plot
figure(4)
[lambda_grid, gain_grid] = meshgrid(lambda_range,gain_range);
mesh(lambda_grid, gain_grid, ssr_matrix') 
xlabel("lambda")
ylabel("gain")
zlabel("ssr")
xlim([lambda_range(1),lambda_range(length(lambda_range))])
ylim([gain_range(1),gain_range(length(gain_range))])
zlim([minimum,max(ssr_matrix(:))])
grid on
title("good model parameter grid")
% fprintf(sprintf("'grid' lambda:"+grid_lambda+"\n'grid' gain:"+grid_gain+"\n"));

y_cap = getGoodModelYCap(u,interval,D,grid_gain,grid_lambda);

ssr = getSsr(deaths,y_cap);
figure(5)
plot(dates, deaths)
title("good model")
hold on
plot(dates, y_cap)
legend("deaths","estimated deaths")

%% functions

function [ssr] = getSsr(y, y_cap)
    ssr = (y-y_cap)'*(y-y_cap);
end

function [y_cap] = getGoodModelYCap(u,interval,D,mu,lambda)
    y_cap = zeros(interval(length(interval)),1);
    for t = interval(1):interval(length(interval))
        for k = 1:t-D-1
            y_cap(t) = y_cap(t) + (u(t-D-k)*mu*lambda*exp(-lambda*k));
        end
    end
    y_cap = y_cap(interval,1);
end