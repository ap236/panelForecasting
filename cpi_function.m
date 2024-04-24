function ppt6cpi (savePath)

% usage: ppt6cpi(savePath)
%
% Calculates forecasts and MSFEs for Pesaran, Pick, Timmermann (2024)
% input is string with path to a folder in which files can be saved

% Andreas Pick

clc;
cc = clock; disp(['Start time: ' num2str(cc(4)) ':' num2str(cc(5))])

% 1) --- Choices that can be made -----

saveOutput = 1;   % 1 - if save output

priorPw = [6, 1, .1]; % prior vector for hierarchical Bayesian forecasts
savePrms = 0; % whether to save all parameters for hierarchical Bayesian forecasts

rollwinVec = 60;  % size of estimation window
setNFore   = 0;   % set number of forecasts, if 0 use max possible

modelVec = [1 4 5]; % choose models to run (1=AR, 2=AR+dfy, 3=AR+tms, 4=AR+PC, 5=AR+all)

burnin = 500; % should be 500
MCMCreps = 1000 + burnin; % should be 1000 + burnin
if savePrms == 0
  savePrms = burnin;
end 

checkEnoughObs = 1;

% --- end of choices ---

% lags in the model:
difflag = 1; % 1 - m-o-m growth date 
lagsCell{1} = [1 2 12]; % lag 1, 2, and 12
nlags = length(lagsCell);

maxposslags      = 12; % maximum possible lags for which to allow presample
maxRollWinLength = 60;

minNoObs      = rollwinVec;  % min no of observations for each individual

nMethods = 10;        % number of methods
benchmark = 3;        % benchmark method: 1 for pooled, 3 for ind.spec.

others = 1:nMethods;
others(benchmark) = [];

% 2) --- load and transform data -----

load('CPIsubindices.mat'); % matrix cpi: subindices in levels
load('CPIsubindices_dates.mat'); % vector dates for subindices in levels
% use data from 1967M1 onwards (as then at least 101 subseries)
cpi = cpi(dates>1967,:);
cpi = cpi(1:end-5,:);      % has 5 obs in 2023 that we don't have macro data for
dates = dates(1:end-5);    

macroPredictors = csvread('MacroPredictors.csv',1,5); % default-spread and term-spread
macroPredictors = macroPredictors(1:912,:); % for some reason Matlab appends rows of 0's that need removing
macroPredictors = macroPredictors(dates>1967,:);

dates = dates(dates>1967);

% calc y-o-y inflation and remove presample from other series

dcpi = 100*(cpi(1+difflag:end,:)./cpi(1:end-difflag,:) - 1); % calculate y-o-y inflation
dfy = macroPredictors(difflag:end-1,1)*100; % already lagged one period compated to CPI
tms = macroPredictors(difflag:end-1,2)*100;
dates = dates(1+difflag:end);

PriceData = dcpi'; 
PriceDates = dates;

% 3) --- forecasting -----

nRegs = 1 + 3 + 1; % AR, 3 regressors, all regressors

for rollwin = rollwinVec % over possibly rolling window sizes

  first_fore_date_location = rollwin+1;
  startLaterThisWin = maxRollWinLength-rollwin;
  if minNoObs > rollwin
    minNoObs = rollwin;
    disp('careful: min number of observations reset to rolling window size')
  end
  factorsCell = cell(size(PriceData,2),1); % reserving space for the PCA factors

  PriceDataFinal = PriceData(:,startLaterThisWin+maxposslags+1:end)';
  dfyFinal = dfy(startLaterThisWin+maxposslags:(end-1));
  tmsFinal = tms(startLaterThisWin+maxposslags:(end-1));
  DatesFinal = PriceDates(startLaterThisWin+maxposslags+1:end);
  [T, N] = size(PriceDataFinal);
  y = PriceDataFinal(:);
  tVec = (1:T)';

  for ilag = 1:nlags 
    lags = lagsCell{ilag};

    y_lagged = nan(size(y,1),length(lags));
    lctr = 0;
    for il = lags
      lctr = lctr+1;
      PriceDataLaggedFinal = PriceData(:,startLaterThisWin+maxposslags-il+1:(end-il))';
      y_lagged(:,lctr) = PriceDataLaggedFinal(:);
    end

    x_vars =  [y_lagged kron(ones(N,1), dfyFinal) kron(ones(N,1), tmsFinal)];

    if setNFore > 0
      nForecasts = setNFore;
    else
      nForecasts = size(DatesFinal(first_fore_date_location:end),1);
    end
    ForecastDates = DatesFinal(first_fore_date_location:end);

    if nForecasts ~= size(DatesFinal(first_fore_date_location:end),1)
      disp(['number of forecasts is set to ' num2str(nForecasts) '!'])
    end


    for iReg = modelVec 

      disp('rollwin ilag iReg')
      disp([rollwin ilag iReg]) % progress info on screen

      if iReg >= (nRegs-1)
        getFactor = 1;
      else
        getFactor = 0;
      end

      % Get rhs variables for this loop (common factor added later):
      if iReg == 1 % pure AR
        X = x_vars(:,1:length(lags));
      elseif iReg < (nRegs-1) % AR with one regressor
        X = x_vars(:,[1:length(lags) length(lags)+iReg-1]);
      elseif iReg == (nRegs-1) % AR with common factors
        X = x_vars(:,1:length(lags));
      elseif iReg == nRegs % % AR with all x
        X = x_vars;
      end

      K = size(X,2)+1+(iReg>3);

      % === priors ====
      priors{1} = zeros(K,1); % prior mean for mean of beta
      priors{2} = eye(K)*(10^(-priorPw(1))); % prior variance for mean of beta, default = eye(K)*(10^-6)
      priors{3} = eye(K)*(10^priorPw(2)); % prior mean for variance of beta, default = eye(K)
      priors{4} = K; % prior d.o.f. for variance of beta
      priors{5} = priorPw(3); % prior mean for error variance, default = 0.1
      priors{6} = 0.1; % prior d.o.f. for error variance

      priorStr = [num2str(priorPw(1)) '_' num2str(priorPw(2)) '_' num2str(priorPw(3))];

      yhatC = cell(nForecasts);
      yerrC = cell(nForecasts);
      ytrue = nan(N,nForecasts);
      nonans = nan(N,nForecasts);
      truncWeights = nan(nForecasts,2);
      NthisVec = nan(nForecasts,1);
      onBoundary = nan(nForecasts,2);
      comps = nan(nForecasts,2,2);

      for tfore = 1:nForecasts 

        yxtvec = (tVec < first_fore_date_location-1+tfore).*(tVec > tfore-1);
        yxtvec = kron(ones(N,1),yxtvec);
        yi = y(yxtvec == 1);
        xi = X(yxtvec == 1,:);

        if getFactor == 1
          % calculate factors
          yFactorin = PriceData(:,1:startLaterThisWin+maxposslags+rollwin)';
          PCAone = principalComponentsNaN(yFactorin,1);
          PCAone = PCAone/sqrt(var(PCAone));
          f_var = PCAone(end-rollwin:end-1);
          xi = [xi kron(ones(N,1),f_var)];
        end

        yxforevec = kron(ones(N,1), (tVec==first_fore_date_location-1+tfore));
        xfore = X(yxforevec==1,:);
        if getFactor == 1
          xfore = [xfore kron(ones(N,1),f_var(end))];
        end
        ytruet = y(yxforevec==1);
        ytrue(:,tfore) = y(yxforevec==1);

        % Get balenced panel by only using units with complete obs in window
        outstrbal = 'Results for balanced panel';
        nany = sum(isnan(reshape(yi,rollwin,N)));
        nanx = zeros(1,N);
        for ik = 1:(K-1)
          nanx = nanx + sum(isnan(reshape(xi(:,ik),rollwin,N)));
        end
        nanyf = isnan(ytrue(:,tfore));
        nanxf = sum(isnan(xfore),2);
        nonans(:,tfore) = (nany+nanx+nanyf'+nanxf' == 0)'; % 1xN vector of 1 if no miss val; else 0
        seli = kron(nonans(:,tfore),ones(rollwin,1));
        ytruet = ytruet(nonans(:,tfore) == 1);

        yibal = yi(seli == 1);
        xibal = xi(seli == 1,:);
        xforebal = xfore(nonans(:,tfore) == 1,:);
        Nthis = length(xforebal);
        Tthis = length(yibal)/Nthis;
        NthisVec(tfore) = Nthis;

        % Panel forecasts ===================
        [yhatF, prmsInd, ~, ~, boundary, components] = panelRandomCorrCoeffFreqForecasts (yibal, xibal, xforebal, Nthis, minNoObs, checkEnoughObs); % unused output: prms2, w, R2
        yhatEB = panelEmpiricalBayesForecasts (yibal, xibal, xforebal, Nthis, minNoObs, checkEnoughObs);
        % starting values
        startVal{1} = prmsInd'; % individual betas
        startVal{2} = mean(prmsInd,2); % mean of betas
        uInd = nan(Nthis*Tthis,1);
        Sig = zeros(K,K);
        for ic = 1:Nthis
          Sig = Sig + (prmsInd(:,ic) - startVal{2})*(prmsInd(:,ic) - startVal{2})'/Nthis;
          uInd((ic-1)*Tthis+1:ic*Tthis) = yibal((ic-1)*Tthis+1:ic*Tthis) - [ones(Tthis,1) xibal((ic-1)*Tthis+1:ic*Tthis,:)]*prmsInd(:,ic); 
        end
        startVal{3} = inv(Sig); % inv_Sigma
        startVal{4} = uInd'*uInd/(Nthis*Tthis-K); % sig2 
        [yhatB, ~] = panelHierarchBayesForecasts (yibal, xibal, xforebal, Nthis, minNoObs, MCMCreps, priors, startVal, savePrms);

        yhatC{tfore} = [yhatF yhatEB yhatB']; 
        yerrC{tfore} = yhatC{tfore} - ytruet*ones(1,nMethods);
        onBoundary(tfore,:) = boundary;
        comps(tfore,:,:) = components;

      end

      % now sort again in order that they originally were
      yhat = nan(N,nForecasts,nMethods);
      yerr = nan(N,nForecasts,nMethods);
      
      for tfore = 1:nForecasts
        for mi = 1:nMethods
          yhat(nonans(:,tfore)==1,tfore,mi) = yhatC{tfore}(:,mi);
          yerr(nonans(:,tfore)==1,tfore,mi) = yerrC{tfore}(:,mi);
        end
      end

      if saveOutput == 1
        saveString = ['inflationForecasts_' num2str(rollwin) '_' num2str(ilag) '_' num2str(iReg)];
        save([savePath saveString],'yerr','yhat','ytrue','nonans','truncWeights','lags','modelVec','maxRollWinLength','lagsCell','MCMCreps','priors','onBoundary','comps');
        disp(saveString)
      end
    end

    cc = clock; disp(['End time calcs: ' num2str(cc(4)) ':' num2str(cc(5))])

  end
end

end % endfunction 
