clear; clc; close all;

figure,
filename = {'T.bin', 'W.bin', 'S.bin'};
plotname = {'Temperature Field T', 'Vorticity W', 'Stream Function S'};

for i = 1 : length(filename)
    fid = fopen(filename{i});
    nx = fread(fid, 1, 'int32');
    ny = fread(fid, 1, 'int32');
    T = reshape(fread(fid, nx*ny, 'double'), nx, ny);
    fclose(fid);
    subplot(length(filename)+1, 1, i), contourf(T'), axis equal;
    title(plotname{i});
end

fid = fopen('v.bin');
t = [];
v = [];
while true
    val = fread(fid, 2, 'double');
    if isempty(val)
        break;
    end
    t(end+1) = val(1);
    v(end+1) = val(2);
end
fclose(fid);
subplot(length(filename)+1, 1, length(filename)+1), plot(t, v), grid;
xlabel('Time'), ylabel('v_{max}');
