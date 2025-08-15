function msfei_i = MCpanelARX (MCrep, N, T, rho0, arho, corrxb, g0, sig2alpha, sig2gamma, MCpart, doSuppl, dataPath)

% Monte Carlo for panel data forecasting project, frequentist forecasts
% called from MCpanelARX_wrapper.m
%
% Input: MCrep - number of repititions for Monte Carlo
%        various parameters

% Andreas Pick - Oct 2024

h = 1; % forecast horizon
Rbx = corrxb*0.5;

checkEnoughObs = 0;

if MCpart == 2
  varyterm = [ones(N/2,1); -ones(N/2,1)];
elseif MCpart == 1
  varyterm = zeros(N,1);
end

% === Set up parameters for Monte Carlo in cell to import to DGP function
pin{1} = rho0;
pin{2} = arho;
pin{3} = g0;
pin{4} = sig2alpha;
pin{5} = sig2gamma;
pin{6} = Rbx;
pin{7} = 1;
pin{8} = 1;
pin{9} = varyterm;

% === Settings to be saved later ===
settingsCell{1,1} = 'MCrep';    settingsCell{1,2} = MCrep;
settingsCell{2,1} = 'N';         settingsCell{2,2} = N;
settingsCell{3,1} = 'T';         settingsCell{3,2} = T;
settingsCell{4,1} = 'rho0';      settingsCell{4,2} = rho0;
settingsCell{5,1} = 'arho';      settingsCell{5,2} = arho;
settingsCell{6,1} = 'corrxb';    settingsCell{6,2} = corrxb;
settingsCell{7,1} = 'g0';        settingsCell{7,2} = g0;
settingsCell{8,1} = 'sig2alpha'; settingsCell{8,2} = sig2alpha;
settingsCell{9,1} = 'sig2gamma'; settingsCell{9,2} = sig2gamma;
settingsCell{10,1} = 'MCpart';   settingsCell{10,2} = MCpart;
settingsCell{11,1} = 'Distribution w_{i,T+1} +/- variance'; settingsCell{11,2} = varyterm;

% === Setting up file name for saving results

savePath = [dataPath 'results/MonteCarlo/'];

if ((g0 > 0) && (rho0>0))
  arx = 1;
  fileStryerr = [savePath 'PPTMCARX_Reps_'];
elseif (rho0>0)
  arx = 0;
  fileStryerr = [savePath 'PPTMCAR_Reps_'];
end

if (rho0 > 0)
  fileStryerr = [fileStryerr ...
    num2str(MCrep) '_N_' num2str(N) ...
    '_T_' num2str(T) '_rho0_' num2str(rho0) '_sig2alpha_' ...
    num2str(sig2alpha) '_cor_beta_x_' ...
    num2str(corrxb)  '_MCpart_' ...
    num2str(MCpart) ];
else
  fileStryerr = [fileStryerr ...
    num2str(MCrep) '_N_' num2str(N) ...
    '_T_' num2str(T) '_g0_' num2str(g0) '_sig2alpha_' ...
    num2str(sig2alpha) '_sig2gam_' ...
    num2str(sig2gamma) '_cor_beta_x_' num2str(corrxb) '_MCpart_' ...
    num2str(MCpart)];
end
if doSuppl == 2
  fileStryerr = [fileStryerr '_wSuppl'];
end
if doSuppl == 1
  fileStryerr = [fileStryerr '_wBayes'];
end


if exist(fileStryerr) ~= 0
  disp([' Already done: ' fileStryerr])
  msfei_i = nan;
  return
end

nMethods_nonSup = 8;
nMethods = nMethods_nonSup + (doSuppl == 2)*4 + (doSuppl == 1)*3; % total number of methods

minNoObs = 1; % not needed here, just for completeness

if g0 > 0
  K = 3;
else
  K = 2;
end

msfei_i = zeros(N,nMethods);
% === Draw parameters ===
seedValueInit = 10000; % seed for rand and randn
parameters = randPanelData16 (N, T, h, pin, seedValueInit, 1);
wTp1 = parameters{10};

if doSuppl == 2
  SupplStuff{1} = parameters{1};
  SupplStuff{2} = parameters{3};
  SupplStuff{3} = parameters{5};
end
if doSuppl == 1
  prior_powers = [6 2 0; 1 2 0];
  priors{1} = zeros(K,1); % prior mean for mean of beta
  priors{2} = eye(K);     % prior precision for mean of beta
  priors{3} = eye(K);     % prior mean for variance of beta
  priors{4} = K;          % prior d.o.f. for variance of beta
  priors{5} = 0.1;        % prior mean for error variance
  priors{6} = 0.1;        % prior d.o.f. for error variance
  MCMCreps = 1500;
end
% === Loop over Monte Carlo replications ===
for imc = 1:MCrep

  % === draw data given parameters
  seedValueInitimc = seedValueInit+imc;
  data = randPanelData16 (N, T, h, parameters, seedValueInitimc, 2);

  y   = vec_ap(data{1}(1:end-h,:));
  if K == 3
    z   = [vec_ap(data{2}(1:end-h,:)) vec_ap(data{3}(1:end-h,:))];
  elseif K == 2
    z   = vec_ap(data{2}(1:end-h,:));
  end
  yf  = data{1}(end,:)'; % Nx1 vector of true values of y to be forecast
  zf = wTp1;

  k = size(z,2);
  if K ~= k+1
    disp(['k ' num2str(k) ' vs K ' num2str(K)])
  end
  K = k+1;

  % === FORECASTS ==================

    [yhat, ~] = panelRandomCorrCoeffFreqForecasts (y, z, zf, N, minNoObs, checkEnoughObs);
    [yhatOG, ~, ~] = panelRandomCoeffForecastsOldWeights (y, z, zf, N, minNoObs, checkEnoughObs);
    yhatEB = panelEmpiricalBayesForecasts (y, z, zf, N, minNoObs, checkEnoughObs);

  if doSuppl == 2
    SupplStuff{4} = data{4};
    yhatSup = panelRandomCorrCoeffFreqForecasts_wOracle (y, z, zf, N, SupplStuff);
  end

  if doSuppl == 1
    % --- hierarchical Bayesian: ---
    prmsInd = nan(K,N);
    for ic = 1:N
      zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
      yi = y((ic-1)*T+1:ic*T);
      % nansi = sum(isnan([zi yi]),2);
      % zi = zi(nansi==0,:);
      % yi = yi(nansi==0);
      prmsInd(:,ic) = pinv(zi'*zi)*(zi'*yi);
    end

    % prmsInd are KxN but bayesian function assumes NxK
    startVal{1} = prmsInd'; % individual betas
    startVal{2} = mean(prmsInd,2); % mean of betas
    uInd = nan(N*T,1);
    Sig = zeros(K,K);
    for ic = 1:N
      Sig = Sig + (prmsInd(:,ic) - startVal{2})*(prmsInd(:,ic) - startVal{2})'/N;
      uInd((ic-1)*T+1:ic*T) = y((ic-1)*T+1:ic*T) - [ones(T,1) z((ic-1)*T+1:ic*T,:)]*prmsInd(:,ic);
    end
    startVal{3} = inv(Sig); % inv_Sigma
    startVal{4} = uInd'*uInd/(N*T-K); % sig2

    for ii = 1:3
      priorsi = priors;
      priorsi{2} = priors{2}*10^(-prior_powers(1,ii));
      priorsi{3} = priors{3}*10^prior_powers(2,ii);
      [yhatBayes{ii}, ~] = panelHierarchBayesForecasts (y, z, zf, N, minNoObs, MCMCreps, priorsi, startVal, 500);
    end

  end

  % Forecast errors:
  yerri = yf*ones(1,nMethods_nonSup) - [yhat(:,1:6) yhatOG(:,1) yhatEB];
  msfei_i(:,1:nMethods_nonSup) = msfei_i(:,1:nMethods_nonSup) + (yerri.^2)/MCrep; % ave over MC
  if doSuppl == 2
    yerri = yf*ones(1,4) - [yhat(:,7:8) yhatSup];
    msfei_i(:,nMethods_nonSup+(1:4)) = msfei_i(:,nMethods_nonSup+(1:4)) + (yerri.^2)/MCrep; % ave over MC
  end
  if doSuppl == 1
    yerri = yf*ones(1,3) - [yhatBayes{1}' yhatBayes{2}' yhatBayes{3}'];
    msfei_i(:,nMethods_nonSup+(1:3)) = msfei_i(:,nMethods_nonSup+(1:3)) + (yerri.^2)/MCrep; % ave over MC
  end

end % === END FORECASTS ==================

% === Calculate MSFE first for each method, then relative to individual forecast ===
msfei = zeros(nMethods,1);
for imsfe = [3 1:2 4:nMethods]
  msfei(imsfe) = mean(msfei_i(:,imsfe))./mean(msfei_i(:,3)); % ave over individuals
end

% === Output ===
save(fileStryerr,'msfei','msfei_i','settingsCell','parameters');

end %function
% --- END ---------------------------------------------------------------------
