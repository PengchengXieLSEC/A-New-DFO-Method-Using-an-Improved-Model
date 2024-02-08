filename = 'problems';
fileID = fopen(filename);
C = textscan(fileID, '%s %f');
fclose(fileID);
b = C{2};
b ./ b
