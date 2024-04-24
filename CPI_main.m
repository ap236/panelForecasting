% Analysing inflation data for panel forecasting project

% Andreas Pick
clear;
clc;

savePath =  'set/path/for/saving/files/';

saveMSFEi = 1;       % 1: save MSFEs to be used for density plots
printMSFEtables = 1; % 1: print MSFE results to screen
printQuantiles = 1;  % 1: print quantiles to screen

% ---------------- CHOICES END ---

cpi_function(savePath); % fun forecasts

% keep the below choices unchanged
rollwinVec = 60;  % size of estimation window
quantiles = [0.01 0.05 0.1 0.5 0.9 0.95 0.99];
nQuantiles = length(quantiles);
lagsCell{1,1} = [1 2 12]; % lag 1, 2, and 12
nlags = length(lagsCell);
sumBest2One = 0;
nFore = 599;
N = 187;


nRegsCombs = 1 + 3 + 1; % AR, 3 regressors, all regressors

yerrComb = cell(5,1);
yhatComb = cell(5,1);
for ii = [1 4 5] % different models used

  if ii == 1

    load([savePath 'inflationForecasts_60_1_1'],'-mat'); 

    yerrComb{ii} = yerr/100;
    yhatComb{ii} = yhat/100;

  elseif ii == 4

    load([savePath 'inflationForecasts_60_1_4'],'-mat'); 

    yerrComb{ii} = yerr/100;
    yhatComb{ii} = yhat/100;

  elseif ii == 5

    load([savePath 'inflationForecasts_60_1_5'],'-mat'); 

    yerrComb{ii} = yerr/100;
    yhatComb{ii} = yhat/100;

  end
  yerr = yerrComb;
  yhat = yhatComb;
end

maxNoModels = size(yerr{1},3);
benchmark = 3;  
subsetMethods = [1:4 (maxNoModels-1):maxNoModels  5 6]; % add 7 8 for equal weights forecasts in appendix 

nMethods = length(subsetMethods);
reOrder = [1 2 4:nMethods]; 
reOrder2 = [3 reOrder];

modelVec = [1 4 5];
nModels = length(modelVec);

% allocate memory
msfe_ave = nan(1,length(modelVec));
msfe_ave_ratio = nan(nMethods,nModels);
beat_benchmark = nan(nMethods,length(modelVec));
fraction_worst = nan(nMethods,length(modelVec));
fraction_best  = nan(nMethods,length(modelVec));
quantiles_msfe = nan(nMethods,length(quantiles),length(modelVec));
dmstat    = nan(nMethods,nModels);
noObs     = nan(nMethods,nModels,N);
indDMstat = cell(nModels,1);
indDMsumStats = cell(nModels,1);

% calculations
mctr = 0;
for model = modelVec
  mctr = mctr+1;

  yerri = yerr{model}(:,:,subsetMethods);
  % --- MSFE---------------------------------------------------------------
  % Calculate the MSFE per unit
  msfe_per_unit = nan(N,nMethods);
  msfe_per_unit(:,:) = nanmean(yerri.^2,2);
  msfe_ave(mctr) = mean(msfe_per_unit(:,benchmark)); % average msfe of benchmark model
  for iMethod = 1:nMethods
    noObs(iMethod,model,:) = sum(~isnan(yerr{model}(:,:,iMethod)),2);
  end

  msfe_ave_ratio(:,mctr) = mean(msfe_per_unit,1)'./(mean(msfe_per_unit(:,benchmark),1)*ones(nMethods,1));

  % how often beat benchmark/how often best/how often worst
  for n = 1:nMethods
    beat_benchmark(n,mctr) = mean(msfe_per_unit(:,n) < msfe_per_unit(:,benchmark));
  end
  fraction_worst(:,mctr) = mean(msfe_per_unit == max(msfe_per_unit,[],2))';
  fraction_best(:,mctr) = mean(msfe_per_unit == min(msfe_per_unit,[],2))';
  % --- Quantile ----------------------------------------------------------
  for m = 1:nMethods
    sorted_ratio_msfe = sort(msfe_per_unit(:,m)./msfe_per_unit(:,benchmark)); % quantile of abs msfe
    quantiles_msfe(m,:,mctr) = sorted_ratio_msfe(round(quantiles*N));
  end

  % --- DM test statistic -------------------------------------------------
  h = 1;
  for cMethod = 1:nMethods
    if cMethod == benchmark
      continue
    end
    zmat = nan(N,nFore);
    zmat(:,:) = yerri(:,:,cMethod).^2 - yerri(:,:,benchmark).^2;
    z = vec(zmat');
    dmstat(cMethod,mctr) = PanelDieboldMarianoPPP_matlab (z, h, N);

     % DM test per unit
    for ind = 1:N
      indDMstat{model}(ind,cMethod) = modifiedDieboldMariano(zmat(ind,:)',1);
    end
    indDMsumStats{model}(:,cMethod) = [ ...
      sum(indDMstat{model}(:,cMethod) < -1.96) ...
      sum((indDMstat{model}(:,cMethod) > -1.96).*(indDMstat{model}(:,cMethod) < 1.96)) ...
      sum(indDMstat{model}(:,cMethod) > 1.96)];
  end
  % --- Saving msfes for graphs -------------------------------------------
  if saveMSFEi == 1
    msfeiout = nan(nMethods,N);
    noObsi = nan(nMethods,N); % not sure I need this
    msfeiout(:,:) = msfe_per_unit'; % ratios are calculated in the plots function
    save('-ascii',[savePath 'msfei_60_' num2str(model) '.txt'],'msfeiout');
    noObsi(:,:) = noObs(:,model,:);
    save('-ascii',[savePath 'noObs_60_' num2str(model) '.txt'],'noObsi');
  end

end


% === Printing stuff =======================================================
absMSFEmultiplicator = 10^5;

matstr = '%1.3f';
for j = 1:nModels-1
  matstr = [matstr ' & %1.3f'];
end
matstr = [matstr ' \n'];

disp('--- CPI ---')

if printMSFEtables == 1

  disp('Abs MSFE of ind spec forecast')
  printStr = sprintf(matstr,msfe_ave*absMSFEmultiplicator);
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
  matstr = [matstr];
  for j = 1:nModels
    for ii = 1:length(reOrder)
      printStr = sprintf(matstr,quantiles_msfe(reOrder(ii),:,j));
      disp(printStr)
    end
    disp('  ')
  end

end
disp('Panel DM test statistics')
disp(dmstat(reOrder,:)')
disp('no of individually significant DM tests stats')
disp('Model 1')
indDMsumStats{1}(:,reOrder)
disp('Model 4')
indDMsumStats{4}(:,reOrder)
disp('Model 5')
indDMsumStats{5}(:,reOrder)

