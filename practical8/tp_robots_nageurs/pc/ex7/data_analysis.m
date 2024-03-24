% Specify the path to your CSV file
filename = 'trajectory-XXX.csv';

% Read the data from the CSV file into a table
data = readtable(filename);

% Display the contents of the table
% disp(data);

% plotting trajectory
plot(data.x, data.y, "Color",[0, 0, 1], "Marker","o")

