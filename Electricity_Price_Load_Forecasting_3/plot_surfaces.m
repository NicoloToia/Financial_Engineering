clear all
close all
clc

% QR-DNN
QR_DNN = ...
readtable('experiments/tasks/EM_price/QR-DNN/recalib_opt_random_1_2/winkler_scores.csv', ReadRowNames=true, Delimiter=";", DecimalSeparator=",")
% N-DNN
N_DNN = ...
readtable('experiments/tasks/EM_price/N-DNN/recalib_opt_random_1_2/winkler_scores.csv', ReadRowNames=true, Delimiter=";", DecimalSeparator=",")
% JSU-DNN
JSU_DNN = ...
readtable('experiments/tasks/EM_price/JSU-DNN/recalib_opt_random_1_2/winkler_scores.csv', ReadRowNames=true, Delimiter=";", DecimalSeparator=",")
% QR-DNN-ARCSINH
QR_DNN_ARCSINH = ...
readtable('experiments/tasks/EM_price/QR-DNN-ARCSINH/recalib_opt_random_1_2/winkler_scores.csv', ReadRowNames=true, Delimiter=";", DecimalSeparator=",")

% plot the 4 surfaces
quantiles = 0.005:0.005:0.05;
hours =  1:24;

figure
% QR-DNN in blue
surf(quantiles, hours, QR_DNN{:,:}, 'FaceColor', 'b')
hold on
% QR-DNN-ARCSINH in yellow
surf(quantiles, hours, QR_DNN_ARCSINH{:,:}, 'FaceColor', 'y')
% fix the view
view(-45, 45)
% add labels
xlabel('Quantile')
ylabel('Hour')
zlabel('Winkler Score')
title('Winkler Score for different models')
legend('QR-DNN', 'QR-DNN-ARCSINH')
hold off

figure
% N-DNN in red
surf(quantiles, hours, N_DNN{:,:}, 'FaceColor', 'r')
hold on
% JSU-DNN in green
surf(quantiles, hours, JSU_DNN{:,:}, 'FaceColor', 'g')

% fix the view
view(45, 45)

% add labels
xlabel('Quantile')
ylabel('Hour')
zlabel('Winkler Score')
title('Winkler Score for different models')
legend('N-DNN', 'JSU-DNN')
hold off
