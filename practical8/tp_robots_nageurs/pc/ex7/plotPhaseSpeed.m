close all
clear all
clc
speeds1 = [0.6012 , 0.3809 , 0.1806];
speeds2 = [0.5702 , 0.4837 , 0.2419];
speeds3 = [0.3747 , 0.2212 , 0.0379];

phases = [0.5, 0.8, 1.2];


% Plot the first subplot
subplot(1,3,1);
plot(phases, speeds1, "Color",[0, 0, 1], "Marker","o")
title('Amplitude: 35°, Frequency: 1Hz');
xlabel('Phase Lag');
ylabel('Speed [m/s]');

% Plot the second subplot
subplot(1,3,2);
plot(phases, speeds2, "Color",[0, 0, 1], "Marker","o")
title('Amplitude: 35°, Frequency: 1.4Hz');
xlabel('Phase Lag');
ylabel('Speed [m/s]');

% Plot the third subplot
subplot(1,3,3);
plot(phases, speeds3, "Color",[0, 0, 1], "Marker","o")
title('Amplitude: 20°, Frequency: 1Hz');
xlabel('Phase Lag');
ylabel('Speed [m/s]');
