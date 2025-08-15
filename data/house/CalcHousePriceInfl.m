% calculate house price inflation from sources
% written for Octave

% Andreas Pick - June 2023

dataPath = 'your/path/here/'

doHouseprices = 0;
doCPI = 1;

if doHouseprices == 1

datain = dlmread('fmhpi_master_file.csv',',',1,0); % reading in the freddie mac data. strings are 0's but we don't need them.
	% columns are: Year, Month, GEO_Type, GEO_Name, GEO_Code, Index_NSA, Index_SA
	% we need Year, Month, GEO_Code, Index_SA (1,2,5,7)
	% the data set includes averages and then the GEO_Code is 0

din = datain(:,[1,2,5,7]);


% (1) Remove averages (GEO-code = 0)
din = din(din(:,3)>0,:);

% (2) Average to quarterly data
dq = []; 

cntr = 1;
while cntr < size(din,1)
	if din(cntr,2) == 1
		if din(cntr+1,2) == 2 && din(cntr+2,2) == 3 && mean(din(cntr:cntr+2,3)) == din(cntr,3)
			dq = [dq; [din(cntr,1) 1 din(cntr,3) mean(din(cntr:cntr+2,4))]];
		elseif mean(din(cntr:cntr+2,3)) == din(cntr,3)
			disp('Something going on here?')
			disp([cntr  din(cntr:cntr+2,3)' din(cntr:cntr+2,1)'])
		endif
	elseif din(cntr,2) == 4
		if din(cntr+1,2) == 5 && din(cntr+2,2) == 6 && mean(din(cntr:cntr+2,3)) == din(cntr,3)
			dq = [dq; [din(cntr,1) 2 din(cntr,3) mean(din(cntr:cntr+2,4))]];
		elseif mean(din(cntr:cntr+2,3)) == din(cntr,3)
			disp('Something going on here?')
			disp([cntr  din(cntr:cntr+2,3)' din(cntr:cntr+2,1)'])
		endif
	elseif din(cntr,2) == 7
		if din(cntr+1,2) == 8 && din(cntr+2,2) == 9 && mean(din(cntr:cntr+2,3)) == din(cntr,3)
			dq = [dq; [din(cntr,1) 3 din(cntr,3) mean(din(cntr:cntr+2,4))]];
		elseif mean(din(cntr:cntr+2,3)) == din(cntr,3)
			disp('Something going on here?')
			disp([cntr  din(cntr:cntr+2,3)' din(cntr:cntr+2,1)'])
		endif
	elseif din(cntr,2) == 10	
		if din(cntr+1,2) == 11 && din(cntr+2,2) == 12 && mean(din(cntr:cntr+2,3)) == din(cntr,3)
			dq = [dq; [din(cntr,1) 4 din(cntr,3) mean(din(cntr:cntr+2,4))]];
		elseif mean(din(cntr:cntr+2,3)) == din(cntr,3)
			disp('Something going on here?')
			disp([cntr din(cntr:cntr+2,3)' din(cntr:cntr+2,1)'])
		endif
	endif
	cntr = cntr+1;
endwhile

% (3) Reshape into form of Cynthia Yang's data

% Cynthia Yang's data:
data = csvread([dataPath 'data/house/data_main_no_header.csv'],1,1);
msa = data(:,1); % MSA code
region = data(:,2); % region
Norig = 377;
Torig = 160;
msam = reshape(msa,Norig,Torig)';    % MSA
regm = reshape(region,Norig,Torig)'; % region

N = 382;
T = 193;
% reshaping, so that data from a quarter are in a row and data from an msa in a column
dqyear  = reshape(dq(:,1),T,N);   % year
dqquart = reshape(dq(:,2),T,N);   % quarter
dqmsa   = reshape(dq(:,3),T,N);   % MSA code
dqhp    = reshape(dq(:,4),T,N);   % house price

% reordering data
dyear = nan(T,Norig);
dquart = nan(T,Norig);
dmsa = nan(T,Norig);
dhp = nan(T,Norig);
cntr = 0;
for ii = 1:377
	colsel = (dqmsa(1,:)==msam(1,ii));
	if sum(colsel) == 1
		cntr++;
		dyear(:,cntr) = dqyear(:,colsel == 1);
		dquart(:,cntr) = dqquart(:,colsel == 1);
		dmsa(:,cntr) = dqmsa(:,colsel == 1);
		dhp(:,cntr) = dqhp(:,colsel == 1);
	elseif sum(colsel) > 1
		disp('more than one col?')
	endif
endfor

save('-mat',[dataPath 'data/house/HousePriceDataJun2023.mat'],'dyear','dquart','dmsa','dhp');


elseif doCPI == 1

	cpi1 = csvread([dataPath data/house/CPIAUCSL.csv'],1,0);
	cpi1 = cpi1(cpi1(:,1) > 1974.9,:); 
	cpi1(end,:) = []; % last month obs is April
	cpi = kron(eye(size(cpi1,1)/3),ones(1,3))*cpi1/3;

	save ('-mat',[dataPath 'data/house/US_CPI_Jun2023.mat'], 'cpi');

endif

