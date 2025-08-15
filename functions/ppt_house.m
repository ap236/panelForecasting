function ppt_house (dataPath, doBayes)

% Analysing housing market data for panel forecasting project

% Andreas Pick - 2025

priorPw{1} =  [6, 1]; % prior powers for variances/precisions set below
priorPw{2} =  [2, 2];
priorPw{3} =  [0, 0];

savePrms = 0;
burnin = 500; % must be greater than 1 if savePrms = 0
MCMCreps = 1000 + burnin;
if savePrms == 0
  savePrms = burnin;
end

checkEnoughObs = 1;
minNoObs = 20; % min no observations per MSA
rollwin = 60; % size of estimation window

wselindc = 1; % min no of adjacent neighbours for MSA to enter forecast eval
lag = 1; % lag length of AR model
h = 1; % forecast horizon

% number of forecast methods:
nMethods = 10 + doBayes*3;

savePath =  [dataPath 'results/house/'];

% ======= LOAD and TRANSFORM DATA
N = 377;
T = 193;
noFore = 128; % number of forecasts

data = csvread([dataPath 'data/house/data_main_no_header.csv'],0,1);
W100 = csvread([dataPath 'data/house/W100.csv']);

data(1,:) = []; % delete header line
region = data(:,2); % region
regm = reshape(region,N,160)'; % region

% load house price data (same format as before) 'dyear' 'dquarter'
% 'dmsa' and 'dhp' (where d stands for data not difference):
load([dataPath 'data/house/HousePriceDataJun2023.mat']);

% loads 'cpi':
load([dataPath 'data/house/US_CPI_Jun2023.mat']);

msam = dmsa;
hopm = dhp;
regm = [regm; ones(33,1)*regm(end,:)];
cpim = kron(ones(1,N),cpi(:,2));
% house price deflation
rhopm = hopm./cpim;
% growth rates for all
dhopm = log(rhopm(5:end,:)./rhopm(1:end-4,:))*100;
dregm = regm(5:end,:);
lhopm = rhopm(5:end,:);
ddyear = dyear(5:end,1);
T = T-4; % adjust sample size for pre-sample lost

% create regional and country averages
dregave = zeros(T,N);
for i = 1:N
  dregave(:,i) = mean(dhopm(:,dregm(1,:)==dregm(1,i)),2);
end
dcountryave = mean(dhopm,2); % Tx1 vector

dhopW = (W100*dhopm')'; % weighting regressors with weighting matrix

Wsel = W100; % weighting matrix to select MSAs to use for forecast evaluation
Wsel = sum(Wsel>0,2);

dhopW(:,sum(W100,2)==0) = []; % removing MSAs without neighbours
dhopm(:,sum(W100,2)==0) = []; % removing MSAs without neighbours
dregm(:,sum(W100,2)==0) = []; % removing MSAs without neighbours
lhopm(:,sum(W100,2)==0) = [];
dregave(:,sum(W100,2)==0) = [];

N = N - sum(sum(W100,2)==0);
Wsel(sum(W100,2)==0) = [];

MSAperReg = nan(8,1); % Number of MSAs per region
for r = 1:8
  MSAperReg(r) = sum(dregm(1,:)==r);
end

% forecasting loops =======================================================

K = 1 + 2*lag + 2; % intercept, lag, spatial lag, regional ave, country ave

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

% setting up result matrices for forecasts of growth rates and levels:
yerr  = cell(noFore,1); %nan(noFore,N,nMethods);
yfore = cell(noFore,1); %nan(noFore,N,nMethods);
wiT   = cell(1,noFore);
wiTp1 = cell(1,noFore);
wthetaiT = cell(1,noFore);
wthetaiTp1 = cell(1,noFore);
for tt = 1:noFore
  wiT{1,tt} = nan(rollwin,N);
  wiTp1{1,tt} = nan(N,1);
  wthetaiT{1,tt} = nan(rollwin,N);
  wthetaiTp1{1,tt} = nan(N,1);
end
ytrue = nan(noFore,N);

Tt = rollwin;

% Loop over forecasting period

parfor ctr = 1:noFore % t = (T-noFore+1-h):(T-h)

  t = ctr + T - noFore - h;
  yforet = nan(N,nMethods);

  y = vec(dhopm(t-rollwin+1:t,:)); % dependent variable

  % Forecasting: ==================================================

  % (1) AR(1) component
  z = nan(rollwin*N,lag);   % regressors in rolling window
  zf = nan(N,lag);          % regressor in forecast period
  zp = dhopm(t-rollwin+1-h:t-h,:);
  zpf = dhopm(t,:);
  z(:,1) = reshape(zp,numel(zp),1); % lagged dep variable
  zf(:,1) = reshape(zpf,numel(zpf),1);
  yf = dhopm(t+h,:); % realisation in forecast period
  yff = yf'*ones(1,nMethods);

  ytrue(ctr,:) = yf;

  % (2) add spatial lag
  zW = nan(size(z,1),size(z,2)+lag); % set up matrix of regressors
  zWf = nan(size(zf,1),size(zf,2)+lag);
  zW(:,1:size(z,2)) = z; % fill with lagged dep variable
  zWf(:,1:size(zf,2)) = zf;
  zW(:,size(z,2)+1) = vec(dhopW(t-rollwin+1-h:t-h,:)); % lagged dep variable
  zWf(:,size(z,2)+1) = vec(dhopW(t,:));

  % (3) add regional and country averages 
  zR = nan(size(zW,1),size(zW,2)+2); % set up matrix of regressors
  zRf = nan(size(zWf,1),size(zWf,2)+2);
  zR(:,1:size(zW,2)) = zW; % fill with lagged dep variable
  zRf(:,1:size(zWf,2)) = zWf;

  zR(:,size(zW,2)+1) = vec(dregave(t-rollwin+1-h:t-h,:)); % lagged dep variable
  zRf(:,size(zW,2)+1) = vec(dregave(t,:));
  zR(:,size(zW,2)+2) = kron(ones(N,1),dcountryave(t-rollwin+1-h:t-h,:)); % lagged dep variable
  zRf(:,size(zW,2)+2) = kron(ones(N,1),dcountryave(t,:));

  [yhat, prmsInd] = panelRandomCorrCoeffFreqForecasts (y, zR, zRf, N, minNoObs, checkEnoughObs);
  yforet(:,1:8) = yhat;

  [yhatOG, ~, ~] = panelRandomCoeffForecastsOldWeights (y, zR, zRf, N, minNoObs, checkEnoughObs);
  yforet(:,9) = yhatOG(:,1);

  yhatEB = panelEmpiricalBayesForecasts (y, zR, zRf, N, minNoObs, checkEnoughObs);
  yforet(:,10) = yhatEB;

  for i = 1:N
    zRi = zR((i-1)*Tt+1:i*Tt,:);
    icov3 = inv(cov(zRi));
    wiT{1,ctr}(:,i) = diag(zRi*icov3*zRi');
    wiTp1{1,ctr}(i) = zRf(i,:)*icov3*zRf(i,:)';
    wthetaiT{1,ctr}(:,i) = [ones(Tt,1) zRi]*prmsInd(:,i);
    wthetaiTp1{1,ctr}(:,i) = [1 zRf(i,:)]*prmsInd(:,i);
  end

  if doBayes == 1
    % starting values
    startVal = cell(4,1);
    startVal{1} = prmsInd'; % individual betas
    startVal{2} = mean(prmsInd,2); % mean of betas
    uInd = nan(N*Tt,1);
    Sig = zeros(K,K);
    for ic = 1:N
      Sig = Sig + (prmsInd(:,ic) - startVal{2})*(prmsInd(:,ic) - startVal{2})'/N;
      uInd((ic-1)*Tt+1:ic*Tt) = y((ic-1)*Tt+1:ic*Tt) - [ones(Tt,1) zR((ic-1)*Tt+1:ic*Tt,:)]*prmsInd(:,ic);
    end
    startVal{3} = inv(Sig); % inv_Sigma
    startVal{4} = uInd'*uInd/(size(uInd,1)-K); % sig2

    [yhatB, ~] = panelHierarchBayesForecasts (y, zR, zRf,   N, minNoObs, MCMCreps, priors(:,1), startVal, savePrms);
    yforet(:,11) = yhatB;
    [yhatB, ~] = panelHierarchBayesForecasts (y, zR, zRf,   N, minNoObs, MCMCreps, priors(:,2), startVal, savePrms);
    yforet(:,12) = yhatB;
    [yhatB, ~] = panelHierarchBayesForecasts (y, zR, zRf,   N, minNoObs, MCMCreps, priors(:,3), startVal, savePrms);
    yforet(:,13) = yhatB;
  end

  yerrt = yforet(:,:) - yff;

  yerr{ctr} = yerrt;
  yfore{ctr} = yforet;

end % end loop over forecast periods (t)

% ------ Output ---------------------

Nthis = sum(Wsel >= wselindc);

yer = nan(noFore,N,nMethods);
yfor = nan(noFore,N,nMethods);
for t = 1:noFore
  yer(t,:,:) = yerr{t}; 
  yfor(t,:,:) = yfore{t};
end
yerr = yer;
yfore = yfor;

save([savePath 'house_' num2str(doBayes)], ...
  'yerr','yfore','ytrue','wiTp1','wiT','wthetaiTp1','wthetaiT', ...
  'Nthis','priors','MCMCreps','Wsel');

end

function v = vec(m)
v = reshape(m,numel(m),1);
end
