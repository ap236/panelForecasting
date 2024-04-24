function housePrice_function(savePath)

  % Analysing housing market data for panel forecasting project
  % includes pre testing forecasts and Maddala's forecasts
  % and loops over forecast windows

  % Andreas Pick

  saveOutput = 1;

  savePrms = 0;
  burnin = 500; % must be greater than 1 if savePrms = 0
  MCMCreps = 1000 + burnin;
  if savePrms == 0
    savePrms = burnin+(savePrms==1);
  end
  priorPw = [6, 1, .1];

  checkEnoughObs = 1;
  startVal = [];

  rollwin = 60; % size of estimation window
  % 125 (125 m dist), 1 (migration), 2 (?), 3 (?)
  qoq = 4; % year on year growth rate for presample
  noFore = 93; % number of forecasts

  wselindc = 1; % min no of adjacent neighbours for MSA to enter forecast eval
  lag = 1; % lag length of AR model
  hvec = 1; % forecast horizon - set up to be direct forecasts for given differencing of data
  doOldData = 0; % 1: use data of Cynthia Yang; else June 2023 data

  % choices end ==========

  seasDumDum = (qoq==1)*1; % dummy if seasonal dummies should be included

  % number of forecast methods:
  nFreqMethods = 8;
  nBayesMethods = 2;
  noMeth = nFreqMethods + nBayesMethods;
  % benchmark = 3;
%  others = 1:noMeth;
%  others(others == benchmark) = [];
  % number of different forecast models:
  noModels = 2; % number of different models used for forecasting
  minNoObs = 20; % min no observations per MSA

  % ======= LOAD and TRANSFORM DATA

  N = 377;
  T = 160;

  data = csvread('data/data_main_no_header.csv',0,1);
  %periods = data(:,1); % time of observation
  data(1,:) = []; % delete header line
  % msa = data(:,1); % MSA code
  region = data(:,2); % region
  cpi = data(:,3); % cpi
  hp = data(:,4); % house prices
  pop = data(:,5); % population
  inc = data(:,6); % personal income

  % data in TxN matrix, where T=160 quarters and N=377 MSAs
  % msam = reshape(msa,N,T)';    % MSA
  hopm = reshape(hp,N,T)';     % house prices
  regm = reshape(region,N,T)'; % region
  cpim = reshape(cpi,N,T)';    % CPI
  popm = reshape(pop,N,T)';    % population
  incm = reshape(inc,N,T)';    % income


  if doOldData == 0
    N = 377;
    T = 193;

    % Data are for 1975Q1 to 2023Q1

    % load new house price data (same format as before) 'dyear' 'dquarter'
    % 'dmsa' and 'dhp' (where d stands for data not difference):
    %load '~/Dropbox (Erasmus Universiteit Rotterdam)/research/panelForecasting/house prices/HousePriceDataJun2023.mat';
    load 'data/HousePriceDataJun2023.mat';

    % loads 'cpi':
    %load '~/Dropbox (Erasmus Universiteit Rotterdam)/research/panelForecasting/house prices/US_CPI_Jun2023.mat';
    load 'data/US_CPI_Jun2023.mat';

    msam = dmsa;
    hopm = dhp;
    regm = [regm; ones(33,1)*regm(end,:)];
    cpim = kron(ones(1,N),cpi(:,2));

    popm = nan(size(cpim)); % not used just for compatability
    incm = nan(size(cpim)); % not used just for compatability

  end

  % house price deflation
  rhopm = hopm./cpim;
  rincm = incm./cpim;

  W100 = csvread('data/W100.csv'); % data from Cythia Yang's Econometrics Reviews paper
  

  % growth rates for all
  dhopm = log(rhopm(5:end,:)./rhopm(1:end-4,:))*100;
  dpopm = log(popm(5:end,:)./popm(1:end-4,:))*100;
  dincm = log(rincm(5:end,:)./rincm(1:end-4,:))*100;
  dregm = regm(5:end,:);
  T = T-qoq; % adjust sample size for pre-sample lost

  % create regional and country averages
  dregave = zeros(T,N);
  for i = 1:N
    dregave(:,i) = mean(dhopm(:,dregm(1,:)==dregm(1,i)),2);
  end
  dcountryave = mean(dhopm,2); % Tx1 vector

  % choice of weighting matrix
  Wthis = W100;
  dhopW = (Wthis*dhopm')'; % weighting regressors with weighting matrix

  Wsel = Wthis; % weighting matrix to select MSAs to use for forecast evaluation
  Wsel = sum(Wsel>0,2);

  dhopW(:,sum(Wthis,2)==0) = []; % removing MSAs without neighbours
  dhopm(:,sum(Wthis,2)==0) = []; % removing MSAs without neighbours
  dpopm(:,sum(Wthis,2)==0) = [];
  dincm(:,sum(Wthis,2)==0) = [];
  dregm(:,sum(Wthis,2)==0) = []; % removing MSAs without neighbours
  dregave(:,sum(Wthis,2)==0) = [];

  N = N - sum(sum(Wthis,2)==0);
  Wsel(sum(Wthis,2)==0) = [];

  MSAperReg = nan(8,1); % Number of MSAs per region
  for r = 1:8
    MSAperReg(r) = sum(dregm(1,:)==r);
  end

  lagW = lag; % lag for spatial AR

  % forecasting loops =======================================================

  K{1} = 1 + length(lag) + 3*seasDumDum; % plain vanilla
  K{2} = K{1} + length(lagW); % add spatial lags
  K{3} = K{2} + 2; % add spatial lags

  % === priors ====
  K1 = K{2};
  K2 = K{3};
  priors{1,1} = zeros(K1,1); % prior mean for mean of beta
  priors{1,2} = zeros(K2,1); % prior mean for mean of beta
  priors{2,1} = eye(K1)*(10^(-priorPw(1))); % prior variance for mean of beta, default: eye(K1)*(10^-6)
  priors{2,2} = eye(K2)*(10^(-priorPw(1))); % prior variance for mean of beta, default: eye(K2)*(10^-6)
  priors{3,1} = eye(K1)*(10^priorPw(2)); % prior mean for variance of beta, default: eye(K1)
  priors{3,2} = eye(K2)*(10^priorPw(2)); % prior mean for variance of beta, default: eye(K2)
  priors{4,1} = K1; % prior d.o.f. for variance of beta
  priors{4,2} = K2; % prior d.o.f. for variance of beta
  priors{5,1} = priorPw(3); % prior mean for error variance, default = 0.1
  priors{5,2} = priorPw(3); % prior mean for error variance, default = 0.1
  priors{6,1} = 0.1; % prior d.o.f. for error variance
  priors{6,2} = 0.1; % prior d.o.f. for error variance

  for h = hvec

    ctr = 0;
    % setting up result matrices for forecasts of growth rates and levels:
    yerr =  cell(noModels,1);
    yfore = cell(noModels,1);
    onBoundaries = cell(noModels,1);
    for jj = 1:noModels
      yerr{jj} = nan(noFore,N,noMeth);
      yfore{jj} = nan(noFore,N,noMeth);
      onBoundaries{jj} = nan(noFore,2);
    end
    truncWeights = nan(noFore,noModels,1);

    % calc the full sample resid variance - assuming one lag only !!!
    if lag ~= 1
      disp('lag not 1: full sample residual variance not calc on correct model')
    end
    yall = nan((T-1)*N,1);
    xall = nan((T-1)*N,7);
    for jj = 1:N
      y = dhopm(2:end,jj);
      x = nan(T-1,7);   % regressors in rolling window
      x(:,1) = dhopm(1:end-1,jj);
      x(:,2) = dhopW(1:end-1,jj);
      x(:,3) = dregave(1:end-1,jj);
      x(:,4) = dcountryave(1:end-1);
      x(:,5) = dincm(1:end-1,jj); % lagged dep variable
      x(:,6) = dpopm(1:end-1,jj); % lagged dep variable
      x(:,7) = ones(T-1,1);
      yall((jj-1)*(T-1)+1:jj*(T-1)) = y;
      xall((jj-1)*(T-1)+1:jj*(T-1),:) = x;
    end
    %nans = isnan(sum([yall xall],2));

    % Loop over forecasting period
    for t = (T-noFore+1-max(hvec)):(T-max(hvec))

      ctr=ctr+1;
      y = vec(dhopm(t-rollwin+1:t,:)); % dependent varible: house price growth in rolling window

      % Forecasting: ==================================================

      % (1) Spatial model =============================================
      z = nan(rollwin*N,lag);   % regressors in rolling window
      zf = nan(N,lag);          % regressor in forecast period
      ct = 0;
      for p = h:(lag+h-1)
        ct=ct+1;
        zp = dhopm(t-rollwin+1-p:t-p,:);
        zpf = dhopm(t-p+h,:);
        z(:,ct) = reshape(zp,numel(zp),1); % lagged dep variable
        zf(:,ct) = reshape(zpf,numel(zpf),1);
      end
      yf = dhopm(t+h,:); % realisation in forecast period
      yff = nan(1,N,noMeth);
      yff(1,:,:)=yf'*ones(1,noMeth);

      zW = nan(size(z,1),size(z,2)+lagW); % set up matrix of regressors
      zWf = nan(size(zf,1),size(zf,2)+lagW);
      zW(:,1:size(z,2)) = z; % fill with lagged dep variable
      zWf(:,1:size(zf,2)) = zf;
      p = h-1;
      for c = size(z,2)+1:size(z,2)+lagW % add with spacial lags
        p=p+1;
        zW(:,c) = vec(dhopW(t-rollwin+1-p:t-p,:)); % lagged dep variable
        zWf(:,c) = vec(dhopW(t-p+h,:));
      end

      [yhatF, prmsInd, ~, ~, boundaries, ~] = panelRandomCorrCoeffFreqForecasts (y, zW, zWf, N, minNoObs, checkEnoughObs);
      yhatEB = panelEmpiricalBayesForecasts (y, zW, zWf, N, minNoObs, checkEnoughObs);
      % starting values
      startVal{1} = prmsInd'; % individual betas
      startVal{2} = mean(prmsInd,2); % mean of betas
      Tt = size(y,1)/N;
      uInd = nan(N*Tt,1);
      Sig = zeros(K{2},K{2});
      for ic = 1:N
        Sig = Sig + (prmsInd(:,ic) - startVal{2})*(prmsInd(:,ic) - startVal{2})'/N;
        uInd((ic-1)*Tt+1:ic*Tt) = y((ic-1)*Tt+1:ic*Tt) - [ones(Tt,1) zW((ic-1)*Tt+1:ic*Tt,:)]*prmsInd(:,ic);
      end
      startVal{3} = inv(Sig); % inv_Sigma
      startVal{4} = uInd'*uInd/(N*Tt-K{2}); % sig2
      [yhatB, ~] = panelHierarchBayesForecasts (y, zW, zWf,   N, minNoObs, MCMCreps, priors(:,1), startVal, savePrms);

      yfore{1}(ctr,:,1:nFreqMethods) = yhatF;
      yfore{1}(ctr,:,nFreqMethods+1) = yhatEB;
      yfore{1}(ctr,:,nFreqMethods+2) = yhatB;

      onBoundaries{1}(ctr,:) = boundaries;

      yerr{1}(ctr,:,:) = yfore{1}(ctr,:,:) - yff;

      % (3) Adding regional and country averages ======================
      zR = nan(size(zW,1),size(zW,2)+2); % set up matrix of regressors
      zRf = nan(size(zWf,1),size(zWf,2)+2);
      zR(:,1:size(zW,2)) = zW; % fill with lagged dep variable
      zRf(:,1:size(zWf,2)) = zWf;
      p = h;
      zR(:,size(zW,2)+1) = vec(dregave(t-rollwin+1-p:t-p,:)); % lagged dep variable
      zRf(:,size(zW,2)+1) = vec(dregave(t,:));
      zR(:,size(zW,2)+2) = kron(ones(N,1),dcountryave(t-rollwin+1-p:t-p,:)); % lagged dep variable
      zRf(:,size(zW,2)+2) = kron(ones(N,1),dcountryave(t,:));

      [yhat, prmsInd, ~, ~, boundaries, ~] = panelRandomCorrCoeffFreqForecasts (y, zR, zRf, N, minNoObs, checkEnoughObs);
      yhatEB = panelEmpiricalBayesForecasts (y, zR, zRf, N, minNoObs, checkEnoughObs);
      % starting values
      startVal{1} = prmsInd'; % individual betas
      startVal{2} = mean(prmsInd,2); % mean of betas
      uInd = nan(N*Tt,1);
      Sig = zeros(K{3},K{3});
      for ic = 1:N
        Sig = Sig + (prmsInd(:,ic) - startVal{2})*(prmsInd(:,ic) - startVal{2})'/N;
        uInd((ic-1)*Tt+1:ic*Tt) = y((ic-1)*Tt+1:ic*Tt) - [ones(Tt,1) zR((ic-1)*Tt+1:ic*Tt,:)]*prmsInd(:,ic);
      end
      startVal{3} = inv(Sig); % inv_Sigma
      startVal{4} = uInd'*uInd/(N*Tt-K{3}); % sig2
      [yhatB, ~] = panelHierarchBayesForecasts (y, zR, zRf,   N, minNoObs, MCMCreps, priors(:,2), startVal, savePrms);

      yfore{2}(ctr,:,1:nFreqMethods) = yhat;
      yfore{2}(ctr,:,nFreqMethods+1) = yhatEB;
      yfore{2}(ctr,:,nFreqMethods+2) = yhatB;

      onBoundaries{2}(ctr,:) = boundaries;

      yerr{2}(ctr,:,:) = yfore{2}(ctr,:,:) - yff;

    end % end loop over forecast periods (t)

    % ------ Output ---------------------

    Nthis = nan(noModels,1);
    for jj = 1:noModels % average over t
      Nthis(jj) = sum(Wsel >= wselindc);
    end

    if saveOutput == 1
      saveString = ['housePriceForecasts.mat'];
      save([savePath saveString], ...
        'yerr','yfore','truncWeights','Wsel','wselindc','Nthis','priors','MCMCreps','onBoundaries');
    end
  end