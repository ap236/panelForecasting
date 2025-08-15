function ppt_cpi (dataPath, doBayes)

% Loads data and forecasts for CPI application 
% for Pesaran, Pick, Timmermann (2025, QE) 

% Andreas Pick - 2025

cc = clock; disp(['Start time: ' num2str(cc(4)) ':' num2str(cc(5))])

% 1) --- Choices that can be made -----

priorPw{1} =  [6, 1]; % prior powers for variances/precisions set below
priorPw{2} =  [2, 2];
priorPw{3} =  [0, 0];

rollwin = 60;  % size of estimation window
setNFore = 0;   % set number of forecasts, if 0 use max possible

% for MCMC:
burnin = 500;
MCMCreps = 1000 + burnin; % replications of MCMC
savePrms = burnin;

checkEnoughObs = 1;

% lags in the model:
difflag = 1; % 1 - m-o-m growth date or 12 y-o-y growth rate
lags = [1 2 12]; % lag 1, 2, and 12
maxlag = max(lags); % largest lag for presample
minNoObs = rollwin;  % min no of observations for each individual

nMethods = 10+doBayes*3;        % number of methods

% 2) --- load and transform data ------------------------------------------

savePath =  [dataPath 'results/cpi/'];

load([dataPath 'data/cpi/CPIsubindices.mat']); % matrix cpi: subindices in levels
load([dataPath 'data/cpi/CPIsubindices_dates.mat']); % vector dates for subindices in levels
% use data from 1967M1 onwards (as then at least 101 subseries)
cpi = cpi(dates>1967,:);
cpi = cpi(1:end-5,:);         % 5 obs in 2023 that we don't have macro
dates = dates(1:end-5);       % data for

macroPredictors = csvread([dataPath 'data/cpi/MacroPredictors.csv'],1,5); % default-spread and term-spread
macroPredictors = macroPredictors(1:912,:); % for some bizzare reason Matlab appends rows of 0's
macroPredictors = macroPredictors(dates>1967,:);

dates = dates(dates>1967);

% calc y-o-y inflation and remove presample from other series

dcpi = 100*(cpi(1+difflag:end,:)./cpi(1:end-difflag,:) - 1); % calculate y-o-y inflation
dfy = macroPredictors(difflag:end-1,1)*100; % already lagged one period compated to CPI
tms = macroPredictors(difflag:end-1,2)*100;
dates = dates(1+difflag:end);

PriceData = dcpi'; % cause that's how the old data was sorted
PriceDates = dates;

% 3) --- forecasting -----

first_fore_date_location = rollwin+1;
startLaterThisWin = 0;
if minNoObs > rollwin
  minNoObs = rollwin;
  disp('careful: min number of observations reset to rolling window size')
end

PriceDataFinal = PriceData(:,startLaterThisWin+maxlag+1:end)';
dfyFinal       = dfy(startLaterThisWin+maxlag:(end-1));
tmsFinal       = tms(startLaterThisWin+maxlag:(end-1));
DatesFinal     = PriceDates(startLaterThisWin+maxlag+1:end);

[T, N] = size(PriceDataFinal);
y = PriceDataFinal(:);
tVec = (1:T)';

y_lagged = nan(size(y,1),length(lags));
lctr = 0;
for il = lags
  lctr = lctr+1;
  PriceDataLaggedFinal = PriceData(:,startLaterThisWin+maxlag-il+1:(end-il))';
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

% Get rhs variables for this loop (common factor added later):
X = x_vars;
K = size(X,2)+2;

% === priors ====
priors{1,1} = zeros(K,1); % prior mean for mean of beta
priors{1,2} = zeros(K,1); % prior mean for mean of beta
priors{1,3} = zeros(K,1); % prior mean for mean of beta
priors{2,1} = eye(K)*(10^(-priorPw{1}(1))); % prior variance for mean of beta, default = eye(K)*(10^-6)
priors{2,2} = eye(K)*(10^(-priorPw{2}(1))); % prior variance for mean of beta, default = eye(K)*(10^-6)
priors{2,3} = eye(K)*(10^(-priorPw{3}(1))); % prior variance for mean of beta, default = eye(K)*(10^-6)
priors{3,1} = eye(K)*(10^priorPw{1}(2)); % prior mean for variance of beta, default = eye(K)
priors{3,2} = eye(K)*(10^priorPw{2}(2)); % prior mean for variance of beta, default = eye(K)
priors{3,3} = eye(K)*(10^priorPw{3}(2)); % prior mean for variance of beta, default = eye(K)
priors{4,1} = K; % prior d.o.f. for variance of beta
priors{4,2} = K; % prior d.o.f. for variance of beta
priors{4,3} = K; % prior d.o.f. for variance of beta
priors{5,1} = 0.1; % prior mean for error variance
priors{5,2} = 0.1; % prior mean for error variance
priors{5,3} = 0.1; % prior mean for error variance
priors{6,1} = 0.1; % prior d.o.f. for error variance
priors{6,2} = 0.1; % prior d.o.f. for error variance
priors{6,3} = 0.1; % prior d.o.f. for error variance

yhatC    = cell(nForecasts,1);
yerrC    = cell(nForecasts,1);
ytrue    = nan(N,nForecasts);
nonansC  = cell(nForecasts,1);

wiT        = cell(nForecasts,1);
wiTp1      = cell(nForecasts,1);
wthetaiT   = cell(nForecasts,1);
wthetaiTp1 = cell(nForecasts,1);

parfor tfore = 1:nForecasts

  yxtvec = (tVec < first_fore_date_location-1+tfore).*(tVec > tfore-1);
  yxtvec = kron(ones(N,1),yxtvec);
  yi = y(yxtvec == 1);
  xi = X(yxtvec == 1,:);

  % calculate factors
  yFactorin = PriceData(:,1:startLaterThisWin+maxlag+rollwin)';
  PCAone = principalComponentsNaN(yFactorin,1);
  PCAone = PCAone/sqrt(var(PCAone));
  f_var = PCAone(end-rollwin:end-1);
  xi = [xi kron(ones(N,1),f_var)];

  yxforevec = kron(ones(N,1), (tVec==first_fore_date_location-1+tfore));
  xfore = X(yxforevec==1,:);
  xfore = [xfore kron(ones(N,1),f_var(end))];
  ytruet = y(yxforevec==1);
  ytrue(:,tfore) = y(yxforevec==1);

  % Get balenced panel by only using units with complete obs in window
  nany = sum(isnan(reshape(yi,rollwin,N)));
  nanx = zeros(1,N);
  for ik = 1:(K-1)
    nanx = nanx + sum(isnan(reshape(xi(:,ik),rollwin,N)));
  end
  nanyf = isnan(ytrue(:,tfore));
  nanxf = sum(isnan(xfore),2);
  nonansi = (nany+nanx+nanyf'+nanxf' == 0)'; % 1xN vector of 1 if no miss val; else 0
  seli = kron(nonansi,ones(rollwin,1));
  ytruet = ytruet(nonansi == 1);

  yibal = yi(seli == 1);
  xibal = xi(seli == 1,:);
  xforebal = xfore(nonansi == 1,:);
  Nthis = length(xforebal);
  Tthis = length(yibal)/Nthis;

  % Panel forecasts ===================
  [yhatF, prmsInd] = panelRandomCorrCoeffFreqForecasts  ...
    (yibal, xibal, xforebal, Nthis, minNoObs, checkEnoughObs);
  [yhatOG, ~, ~] = panelRandomCoeffForecastsOldWeights (yibal, xibal, ...
    xforebal, Nthis, minNoObs, checkEnoughObs);
  yhatEB = panelEmpiricalBayesForecasts (yibal, xibal, xforebal, Nthis, ...
    minNoObs, checkEnoughObs);
  yhatC{tfore,1} = [yhatF yhatOG(:,1) yhatEB];
  if doBayes == 1
    rng( tfore, 'twister');
    % starting values
    startVal = cell(4,1);
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
    % prior 1
    [yhatB, ~] = panelHierarchBayesForecasts (yibal, xibal, xforebal, Nthis, minNoObs, MCMCreps, priors(:,1), startVal, savePrms);
    yhatC{tfore,1} = [yhatC{tfore,1} yhatB'];
    % prior 2
    [yhatB, ~] = panelHierarchBayesForecasts (yibal, xibal, xforebal, Nthis, minNoObs, MCMCreps, priors(:,2), startVal, savePrms);
    yhatC{tfore,1} = [yhatC{tfore,1} yhatB'];
    % prior 3
    [yhatB, ~] = panelHierarchBayesForecasts (yibal, xibal, xforebal, Nthis, minNoObs, MCMCreps, priors(:,3), startVal, savePrms);
    yhatC{tfore,1} = [yhatC{tfore,1} yhatB'];
  end
  yerrC{tfore,1} = yhatC{tfore,1} - ytruet*ones(1,nMethods);
  nonansC{tfore,1} = nonansi;

  wiT{tfore} = nan(rollwin,N);
  wiTp1{tfore} = nan(N,1);
  wthetaiT{tfore} = nan(rollwin,N);
  wthetaiTp1{tfore} = nan(N,1);
  ictr = 0;
  for ii = 1:N
    if nonansi(ii) == 1
      ictr = ictr+1;
      xii = xibal((ictr-1)*Tthis+1:ictr*Tthis,:);
      icov = inv(cov(xibal));
      wiT{tfore}(:,ii) = diag(xii*icov*xii');
      wiTp1{tfore}(ii) = xforebal(ictr,:)*icov*xforebal(ictr,:)';
      wthetaiT{tfore}(:,ii) = [ones(Tthis,1) xii]*prmsInd(:,ictr);
      wthetaiTp1{tfore}(ii) = [1 xforebal(ictr,:)]*prmsInd(:,ictr);
    end
  end
end

% now sort again in order that they originally were
yhat = nan(N,nForecasts,nMethods);
yerr = nan(N,nForecasts,nMethods);

for tfore = 1:nForecasts
  for mi = 1:nMethods
    yhat(nonansC{tfore,1}==1,tfore,mi) = yhatC{tfore}(:,mi);
    yerr(nonansC{tfore,1}==1,tfore,mi) = yerrC{tfore}(:,mi);
  end
end

saveString = ['cpi_forecasts_doBayes_' num2str(doBayes)];
save([savePath saveString],'yerr','yhat','ytrue','wiT','wiTp1','wthetaiT','wthetaiTp1','nonansC','lags','MCMCreps','priors');