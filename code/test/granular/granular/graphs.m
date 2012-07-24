
inputfile = 'dump';


output_file_name = [inputfile,'.dat'];


fid = fopen(output_file_name);
data = fread(fid, inf, 'float64');
fclose(fid);

body_count = 1000;
channels = 3 * body_count;

timesteps = floor(length(data)/channels);
%throw away incomplete timesteps at the end
data = data(1:(timesteps * channels));
trajectory = reshape(data, channels, timesteps);

i = 1:timesteps;

x0 = trajectory(1,:);
y0 = trajectory(2,:);
z0 = trajectory(3,:);
x1 = trajectory(4,:);
y1 = trajectory(5,:);
z1 = trajectory(6,:);
x2 = trajectory(7,:);
y2 = trajectory(8,:);
z2 = trajectory(9,:);
zi = trajectory(999 * 3,:);
%...
    
plot(i, zi);