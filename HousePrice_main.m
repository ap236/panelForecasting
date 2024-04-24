% Analysing housing market data for panel forecasting project
% simplified

% Andreas Pick 
clear;
clc;

saveMSFEi = 1;          % save MSFEi for further analysis (density plots)
printMSFEtables = 1;    % print MSFE results on screen
printQuantiles = 1;     % print quantiles on screen

savePath =  'path/to/save/files/to'; % change the path to one that works for you

housePrice_function(savePath);

% some choices ========

doFreq = 1;
doBayes = 1;

quantiles = [0.01 0.05 0.10 0.50 0.90 0.95 0.99];
nQuantiles = length(quantiles);
quantilesNarrow = [0.05 0.10 0.50 0.90 0.95];
nQuantilesNarrow = length(quantilesNarrow);

% choices end ==========

% settings that don't change:
nFore = 93; % number of forecasts per MSA
qoq = 4; % if 1: quarter on quarter growth, otherwise year on year
seasDumDum = 0;
sumBest2One = 0;  % ensure that best forecasts methods sum to one
% (if several same, they get proportional score)
medianMSFE = 0; % if 1 report the median of the individual MSFEs
wselindc = 1; % min no of adjacent neighbours for MSA to enter forecast eval
lag = 1;      % lag length of AR model
h = 1;        % forecast horizon - set up to be direct forecasts for given
rollwin = 60;

nFreq = 6;               % standard and shrinkage forecasts
nBayes = 2;              % number of Bayesian forecasts

% number of different forecast models:
minNoObs = 20; % min no observations per MSA
nModels = 2; % (1) SAR (2) SARX

N = 377;

% ------ Load Forecasts ---------------------
               
load([savePath 'housePriceForecasts'],'-mat');

maxNoModels = size(yerr{1},3);
benchmark = 3;
subsetMethods = [1:4 (maxNoModels-1):maxNoModels 5 6]; % 7 8 is eq. weights for appendix
others = subsetMethods;
others(others==benchmark)=[];
nMethods = length(subsetMethods);
reOrder = [1 2 4:nMethods];
reOrder2 = [benchmark reOrder];

% ------ Output ---------------------

msfe_ave_benchmark = nan(1,nModels);
msfe_ave = nan(nMethods,nModels);
msfe_per_unit_ratio = nan(nMethods,nModels);
msfe_ave_ratio = nan(nMethods,nModels);
beat_benchmark = nan(nMethods,nModels);
fraction_worst = nan(nMethods,nModels);
fraction_best = nan(nMethods,nModels);
quantiles_msfe = nan(nMethods,nQuantiles,nModels);
quantiles_msfe_per_period = nan(nMethods,nQuantiles,nModels);

dmstat = nan(nMethods-1,nModels);
noObs = nan(nMethods,nModels,N);
indDMstat = cell(nModels,1);
indDMsumStats = cell(nModels,1);

Nthis = nan(nModels,1);
mctr = 0;
for cModel = [1 2] % average over t
  mctr = mctr+1;
  jjy = cModel;
  yerri = yerr{jjy}(:,:,subsetMethods);
  Nthis(cModel) = sum(Wsel >= wselindc);
  nForecasts = size(yerr{jjy},1);
  msfe_per_unit = nan(Nthis(cModel),nMethods);
  msfe_per_period = nan(nForecasts,nMethods);
  for cMethod = [benchmark 1 2 4:nMethods]
    msfe_per_unit(Wsel >= wselindc,cMethod) = nanmean(yerri(:,Wsel >= wselindc,cMethod).^2,1);
    msfe_per_period(:,cMethod) = nanmean(yerr{jjy}(:,Wsel >= wselindc,cMethod).^2,2);
    noObs(cMethod,cModel,Wsel >= wselindc) = sum(~isnan(yerr{jjy}(:,Wsel >= wselindc,cMethod)),1);
    msfe_per_unit_ratio(cMethod,cModel) = nanmean(msfe_per_unit(:,cMethod)./(msfe_per_unit(:,benchmark)),1);
  end
  msfe_ave(:,cModel) = nanmean(msfe_per_unit,1);
  msfe_ave_benchmark(cModel) = nanmean(msfe_per_unit(:,benchmark));

  msfe_per_unit_ratios = msfe_per_unit./(msfe_per_unit(:,benchmark)*ones(1,nMethods));
  msfe_ave_ratio(:,cModel) = mean(msfe_per_unit,1)'./(mean(msfe_per_unit(:,3),1)*ones(nMethods,1));

  % --- DM test statistic -------------------------------------------------
  h = 1;
  indDMstat{cModel} = nan(Nthis(cModel),nMethods);
  indDMsumStats{cModel} = nan(3,nModels);
  dctr = 0;
  for cMethod = subsetMethods
    if cMethod == benchmark
      continue
    end
    dctr=dctr+1;
    % panel DM test of Pesaran and Pick (2013)
    zmat = nan(size(yerr{jjy},1),Nthis(cModel));
    zmat(:,:) = yerr{jjy}(:,Wsel >= wselindc,cMethod).^2 - yerr{jjy}(:,Wsel >= wselindc,benchmark).^2;
    z = vec(zmat);
    dmstat(dctr,cModel) = PanelDieboldMarianoPPP_matlab (z, h, Nthis(cModel));
    % DM test per unit
    for ind = 1:Nthis(cModel)
      indDMstat{cModel}(ind,dctr) = modifiedDieboldMariano(zmat(:,ind),1);
    end
    indDMsumStats{cModel}(:,dctr) = [ ...
      sum(indDMstat{cModel}(:,dctr) < -1.96) ...
      sum((indDMstat{cModel}(:,dctr) > -1.96).*(indDMstat{cModel}(:,dctr) < 1.96)) ...
      sum(indDMstat{cModel}(:,dctr) > 1.96)];
  end


  % how often beat benchmark/how often best/how often worst
  for n = 1:nMethods
    beat_benchmark(n,mctr) = nanmean(msfe_per_unit(:,n) < msfe_per_unit(:,benchmark));
  end
  fraction_worst(:,mctr) = nanmean(msfe_per_unit == max(msfe_per_unit,[],2))';
  fraction_best(:,mctr) = nanmean(msfe_per_unit == min(msfe_per_unit,[],2))';
  % --- Quantile ----------------------------------------------------------
  for m = 1:nMethods
    sorted_ratio_msfe = sort(msfe_per_unit(:,m)./msfe_per_unit(:,benchmark)); % quantile of abs msfe
    sorted_ratio_msfe = sorted_ratio_msfe(~isnan(sorted_ratio_msfe));
    quantiles_msfe(m,:,mctr) = sorted_ratio_msfe(round(quantiles*Nthis(cModel)));
    % and now period quantiles
    sorted_ratio_msfe = sort(msfe_per_period(:,m)./msfe_per_period(:,benchmark)); % quantile of abs msfe
    sorted_ratio_msfe = sorted_ratio_msfe(~isnan(sorted_ratio_msfe));
    quantiles_msfe_per_period(m,:,mctr) = sorted_ratio_msfe(round(quantiles*nForecasts));
  end


  if saveMSFEi == 1
    msfeiout = nan(nMethods,Nthis(cModel));
    msfeiout(:,:) = msfe_per_unit';
    save('-ascii',[savePath 'housePrice_msfei_model_' num2str(cModel) '.txt'],'msfeiout');
  end

end



%% printing to screen ===================================================
absMSFEmultiplicator = 1;

matstr = '%1.3f';
for j = 1:nModels-1
  matstr = [matstr ' & %1.3f'];
end
matstr = [matstr ' \n'];

disp('--- house prices ---')

if printMSFEtables == 1

  disp('Abs MSFE of ind spec forecast')
  printStr = sprintf(matstr,msfe_ave(3,:)*absMSFEmultiplicator);
  disp(printStr)
  disp('Rel MSFE to ind forecast')
  printStr = sprintf(matstr,msfe_ave_ratio(reOrder,:)');
  disp(printStr)

  disp('prop better than ind forecast')
  printStr = sprintf(matstr,beat_benchmark(reOrder,:)');
  disp(printStr)

  disp('Prop best forecast')
  printStr = sprintf(matstr,fraction_best(reOrder2,:)');
  disp(printStr)

  disp('Prop worst forecast')
  printStr = sprintf(matstr,fraction_worst(reOrder2,:)');
  disp(printStr)
end
if printQuantiles == 1
  disp('Quantiles')
  matstr = '%1.3f';
  for j = 1:length(quantiles)-1
    matstr = [matstr ' & %1.3f'];
  end
  for j = 1:nModels
    for ii = 1:length(reOrder)
      printStr = sprintf(matstr,quantiles_msfe(reOrder(ii),:,j));
      disp(printStr)
    end
    disp('  ')
  end
end

disp('Panel DM test statistics')
disp(dmstat(:,1:2)')
disp('no of individually significant DM tests stats')
disp('SAR')
disp(indDMsumStats{1})
disp('SARX')
disp(indDMsumStats{2})
