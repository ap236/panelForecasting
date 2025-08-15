function [yhat, prms] = panelHierarchBayesForecasts (y, x, xf, N, minNoObs, MCMCreps, prior, startVal, savePrms)

  % usage: [yhat, prms] = panelHierarchBayesForecasts (y, x, x_forecast_period, N, minNoObs, MCMCreps, priors, startVal)
  %
  % forecasts using Bayesian hierarchical model
  %
  % y - T*Nx1
  % x - T*NxK, regressors without intercept
  % x_in_forecast_period - NxK, regressors in T+h without intercept
  % N - scaler, number of individuals
  % minNoObs - minimum number of observations for any individual, scalar
  % MCMCreps - number of MCMC iterations
  % priors - cell with priors
  % startVals (optional) - starting values for MCMC
  %
  % yhat - Nx1, forecasts from different forecast methods
  % prms - Nxk parameter matrix

  % Andreas Pick - 2025

  [NT,K] = size(x);
  T = NT/N;
  K = K+1;

  % eliminating individuals without sufficient observations
  enoughObs = zeros(N,1);
  bols = nan(N,K);
  s2ols = nan(N,1);
  siXX = zeros(K,K);
  for ic = 1:N
    xi = [ones(T,1) x((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([xi yi]),2);
    xi = xi(nansi==0,:);
    if size(xi,1) < minNoObs
      continue
    elseif sum(isnan(xf(ic,:))) > 0
      continue
    elseif rank(xi'*xi) < K
      continue
    end
    enoughObs(ic) = 1;
    bols(ic,:) = (xi\yi)';
    u = yi - xi*bols(ic,:)';
    s2ols(ic) = u'*u/(T-K);
    siXX = siXX + s2ols(ic)*inv(xi'*xi);
  end
  x = x(kron(enoughObs,ones(T,1))==1,:);
  y = y(kron(enoughObs,ones(T,1))==1);
  xf = xf(enoughObs==1,:);
  N = sum(enoughObs);
  findObs = find(enoughObs==1);

  if N < 2
    prms = nan;
    return;
  end

  priormb = prior{1}; % prior mean for mean of beta
  prioriSb = prior{2}; % prior precision for mean of beta
  priorVb = prior{3}; % prior mean for variance of beta
  priornb = prior{4}; % prior d.o.f. for variance of beta
  priors2 = prior{5}; % prior mean for error variance
  priorns = prior{6}; % prior d.o.f. for error variance

  if savePrms == 1 % === STORE ALL DRAWS, ALSO FOR PRMS ===
    
    % setting up storage for what I draw
    yhat = nan(N,MCMCreps);
    betas = nan(N,K,MCMCreps);
    mub = nan(K,MCMCreps);
    iSig = nan(K,K,MCMCreps);
    sig2 = nan(MCMCreps,1);
    yhatB = nan(N,MCMCreps);

    if nargin == 7
    	betas(:,:,1) = zeros(N,K); % individual betas
    	mub(:,1) = zeros(K,1); % mean of betas
    	iSig(:,:,1) = eye(k);
    	sig2(1) = 0.1;
    else
    	betas(:,:,1) = startVal{1}; % individual betas
    	mub(:,1) = startVal{2}; % mean of betas
    	iSig(:,:,1) = startVal{3};
    	sig2(1) = startVal{4};
    end

    for rep = 2:MCMCreps

      sumsqres = 0;

    	for ic = 1:N
    		% draw individual betas
    		xi = [ones(T,1) x((ic-1)*T+1:ic*T,:)];
    		yi = y((ic-1)*T+1:ic*T,:);
    		Di = inv(xi'*xi/sig2(rep-1) + iSig(:,:,rep-1));
    		betas(ic,:,rep) = (Di*(xi'*yi/sig2(rep-1) + iSig(:,:,rep-1)*mub(:,rep-1)) + chol(Di)'*randn(K,1))';
    		
    		% calculate forecast
  	  	yhatB(ic,rep) = betas(ic,:,rep) * [1, xf(ic,:)]';  

        sumsqres = sumsqres + sum((yi - xi*betas(ic,:,rep)').^2);

    	end

   	  % draw standard deviations: Gelfand et al.(1990, JASA)
 		  df = (T*N + priorns)/2;
    	s = 1/2*(sumsqres + priorns*priors2);
 		  isig2 = gamrnd(df, 1/s);
      sig2(rep) = 1/isig2;

    	% draw mu_beta
    	barbeta = mean(betas(:,:,rep),1)';
    	V = inv(N*iSig(:,:,rep-1) + prioriSb);
    	mub(:,rep) = V*(N*iSig(:,:,rep-1)*barbeta + prioriSb*priormb) + chol(V)'*randn(K,1);

    	% draw V_beta
    	df = (N+priornb);
      W = (betas(:,:,rep)' - mub(:,rep)*ones(1,N))*(betas(:,:,rep) - ones(N,1)*mub(:,rep)') + priorVb*priornb;
  	  iSig(:,:,rep) = wishrnd(inv(W),df);
    end

    yhat{1}(findObs,:) = yhatB;

    prms{1} = betas;
    prms{2} = mub;
    prms{3} = iSig;
    prms{4} = sig2;

  else % ==== ONLY SAVE FORECASTS, NOT PARAMETERS =============

    % setting up storage for what I draw
    yhatB = zeros(N,1);

    if nargin == 7
    	betas = zeros(N,K); % individual betas
    	mub = zeros(K,1); % mean of betas
    	iSig = eye(k);
    	sig2 = 0.1;
    else
    	betas = startVal{1}; % individual betas
    	mub = startVal{2}; % mean of betas
    	iSig = startVal{3};
    	sig2 = startVal{4};
    end

    for rep = 2:MCMCreps

      sumsqres = 0;
    	% draw beta_i and the sigma^2_i
    	for ic = 1:N
    		xi = [ones(T,1) x((ic-1)*T+1:ic*T,:)];
    		yi = y((ic-1)*T+1:ic*T,:);
    		% draw individual betas
    		Di = inv((xi'*xi)/sig2 + iSig);

    		betas(ic,:) = (Di*(xi'*yi/sig2 + iSig*mub) + chol(Di)'*randn(K,1))';

    		% calculate forecast
        if rep > savePrms
    	  	yhatB(ic) = yhatB(ic) + (betas(ic,:) * [1, xf(ic,:)]')/(MCMCreps-savePrms);  
        end

        sumsqres = sumsqres + sum((yi - xi*betas(ic,:)').^2);

    	end

   	  % draw standard deviations: Gelfand et al.(1990, JASA)
 		  df = (T*N + priorns)/2;
    	s = 1/2*(sumsqres + priorns*priors2);
 		  isig2 = gamrnd(df, 1/s);
      sig2 = 1/isig2;

    	% draw mu_beta
    	barbeta = mean(betas,1)';
    	V = inv(N*iSig + prioriSb);
    	mub = V*(N*iSig*barbeta + prioriSb*priormb) + chol(V)'*randn(K,1);

    	% draw V_beta
    	df = (N+priornb);
      W = (betas' - mub*ones(1,N))*(betas - ones(N,1)*mub') + priorVb*priornb;

  	  iSig = wishrnd(inv(W),df); 
    end

    yhat(findObs) = yhatB;
    prms = [];
  end
end