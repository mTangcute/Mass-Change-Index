% Description:
%   This script computes the GEV-based Mass Change Index (MCI) from
%   line-of-sight gravity difference (LGD) anomalies.
%
%   The LGD anomaly is defined as:
%
%       LGD_anomaly = LGD_measurements - LGD_climatology
%
%   A Generalized Extreme Value (GEV) distribution is fitted to the LGD
%   anomaly series. The fitted cumulative distribution function (CDF) values
%   are then transformed into standard normal quantiles to obtain the
%   GEV-based MCI.
%
% Input:
%   LGD_measurements: LGD measurement time series
%   LGD_climatology:  climatological LGD time series
%   dates:            the corresponding time
%
%   The two time series are calculated by integrating LGD data
%   crossing the target block.
%
% Output:
%   MCI:              Derived Mass Change Index time series
%
% Author: Miao Tang
% Date: 20/04/2026
% Institution: Southwest Jiaotong University
% E-mail: miaotang@my.swjtu.edu.cn

clc;
clear;
close all;

% Read and load the data
input_file = 'mssp_lgd2mci.mat';

data = load(input_file);

LGD_measurements = data.LGD_measurements;
LGD_measurements = LGD_measurements(:) * 1e9; % m/s^-2 -> nm/s^-2

LGD_climatology = data.LGD_climatology;
LGD_climatology = LGD_climatology(:) * 1e9;   % m/s^-2 -> nm/s^-2

dates = data.dates;
dates = dates(:);

% Prepare the LGD anomaly time series
LGD_anomaly_all = LGD_measurements - LGD_climatology;

valid_idx = ~isnan(LGD_anomaly_all) & ...
            ~isinf(LGD_anomaly_all);

LGD_anomaly = LGD_anomaly_all(valid_idx);
dates_valid = dates(valid_idx);

if length(LGD_anomaly) < 20
    error('Too few valid samples for reliable GEV fitting.');
end

% Fit the LGD anomaly time series using the GEV distribution
parmhat = gevfit(LGD_anomaly);

shape_parameter    = parmhat(1);
scale_parameter    = parmhat(2);
location_parameter = parmhat(3);

% Calculate the GEV-based CDF
cdf_values = gevcdf(LGD_anomaly, ...
                    shape_parameter, ...
                    scale_parameter, ...
                    location_parameter);

cdf_values(cdf_values <= 0) = eps;       % Avoid Inf values from norminv(0)
cdf_values(cdf_values >= 1) = 1 - eps;   % Avoid Inf values from norminv(1)

% Calculate the MCI
MCI = norminv(cdf_values);

% Make the figure
figure;
plot(dates_valid, MCI, 'b', 'LineWidth', 1.5);
grid on;
box on;
xlabel('Time');
ylabel('Index');
title('Mass Change Index');