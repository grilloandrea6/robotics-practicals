close all
clear all
clc
speeds1 = [0.6012 , 0.3809 , 0.1806];
err1 = [0.0419, 0.002, 0.0147];
speeds2 = [0.5702 , 0.4837 , 0.2419];
err2 = [0.0051 0.0271 0.0049];
speeds3 = [0.3747 , 0.2212 , 0.0379];
err3 = [0.0221, 0.0318, 0.0442];
phases = [0.5, 0.8, 1.2];


% Plot the first subplot
subplot(1,3,1);
errorbar(phases, speeds1, err1, "Color",[0, 0, 1], "Marker","o")
xlim([0.4 1.3])
title('Amplitude: 35°, Frequency: 1Hz');
xlabel('Phase Lag');
ylabel('Speed [m/s]');

% Plot the second subplot
subplot(1,3,2);
errorbar(phases, speeds2, err2, "Color",[0, 0, 1], "Marker","o")
xlim([0.4 1.3])
title('Amplitude: 35°, Frequency: 1.4Hz');
xlabel('Phase Lag');
ylabel('Speed [m/s]');

% Plot the third subplot
subplot(1,3,3);
errorbar(phases, speeds3,err3, "Color",[0, 0, 1], "Marker","o")
xlim([0.4 1.3])
title('Amplitude: 20°, Frequency: 1Hz');
xlabel('Phase Lag');
ylabel('Speed [m/s]');
