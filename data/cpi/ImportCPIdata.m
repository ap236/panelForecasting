% import the csv files with the CPI data and put them into shape
% Gnu-Octave script
clear
clc

dataPath = "your/data/path";

cpiDirectory = [dataPath "data/cpi/series/"];
fileNames = dir(cpiDirectory);

cpidata = nan(917,188);
names = cell(187,1);

for ci = 3:189 % first two are . and .. so from element 3 to 189

	toread = [cpiDirectory fileNames(ci).name];
	names{ci-2} = substr(fileNames(ci).name,1,length(fileNames(ci).name)-4); 
	d = dlmread(toread,',',1,1);
	years = real(d(:,1));
	data = d(:,3);

	if ci == 3
		cpidata(:,1) = years;
	endif

	missingAtEnd = 5*2023 - sum(years(end-4:end,1)); % there should be 5 month in 2023
	missingAtBeg = 917 - size(data,1) - missingAtEnd; % 

	cpidata(missingAtBeg+1:917-missingAtEnd,ci-1) = data;

end

% put this whole thing into a new csv file
cpiDirectory2 = [dataPath "data/cpi/"];
outfile = fopen([cpiDirectory2 "CPIsubindices.csv"],'a');

dataString = "%4.0f,";
namesString = "Dates,";
for i = 1:186
	namesString = [namesString names{i} ","];
	dataString = [dataString "%4.4f,"];
endfor
namesString = [namesString names{187}];
dataString = [dataString "%4.4f\n"];

fprintf(outfile,"%s\n",namesString);
for i = 1:917
	fprintf(outfile,dataString,cpidata(i,:));
endfor


fclose(outfile);

