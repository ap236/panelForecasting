% Main script
clear; close('all'); clc;

% Choices and things to set prior to running the code:
% Set path
dataPath = '~/your/path/to/PPTPanelForecasting/';

% 1. Applications
doHouse      = 1; % 1: run/tabulate the house price application
doCPI        = 1; % 1: run/tabulate the CPI application
% Do forecasts or only tabulate existing forecast results?
doForecastsHouse = 2; % 1: run forecasts, 2: use previously obtained results
doForecastsCPI   = 2; % 1: run forecasts, 2: use previously obtained results
% For applications: results in text (0) or in supplement (1)
doSupplement = 1; % 0: do not run hierarchical Bayesian model;
                  % 1: run hierarchical Bayesian models

% 2. Monte Carlo
doMonteCarlo = 0; % 1: produce results of the Monte Carlo experiments, 2: only tabulare results from previous runs
% For Monte Carlo: tables in main text (0) or in supplement and, if so,
doMCSupplement = 1; % 1: hierarchical Bayesian forecasts in Supplement
                    % 2: Oracle weights and equal weights in Supplement

% 3. Additional tables
doTableA1 = 0; % 1: produce Table A.1
doTableS1 = 0; % 1: produce Table S.1




% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% Do not alter code below
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% 1. house price application
if doHouse == 1
  if doForecastsHouse == 1
    ppt_house (dataPath, doSupplement); % runs forecasts
  end
  ppt_house_tabulate (dataPath, doSupplement); % analyses forecasts
  if doSupplement == 0
    ppt_house_density_plots (dataPath); % plots forecast results
  end
end
% 2. CPI application
if doCPI == 1
  if doForecastsCPI == 1
    ppt_cpi (dataPath, doSupplement); % runs forecasts
  end
  ppt_cpi_tabulate (dataPath, doSupplement); % analyses forecasts
  if doSupplement == 0
    ppt_cpi_density_plots (dataPath, doHouse); % plots forecast results
  end
end

% Monte Carlo results
if doMonteCarlo > 0
  if doMonteCarlo == 1
    parfor i = 1:72
      if i >= 37
        wtp1 = 2;
      else
        wtp1 = 1;
      end
      if wtp1 == 1
        runno = i;
      else
        runno = i-36;
      end
      MCpanelARX_wrapper(runno, wtp1, doMCSupplement, dataPath);
    end
  end
  MCpanelARX_tabulate (dataPath, doMCSupplement);
end

% Table A.1
if doTableA1 == 1
  disp('--- Table A.1 ---')
  beta0 = [0.3 0.45 0.49 0.4999];
  outp = nan(length(beta0),2);
  for ii = 1:length(beta0)
    [outp(ii,1), outp(ii,2)] = DeltaARnumerical (beta0(ii));
    printString = [num2str(beta0(ii)) ' & %1.3f & %1.3f'];
    disp(sprintf(printString, outp(ii,:)))
  end
end

% Table S.1
if doTableS1 == 1
  getPanelR2Fun;
end