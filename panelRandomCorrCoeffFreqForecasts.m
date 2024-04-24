function [yhat, betaInd, weights, R2, boundaries, components] = panelRandomCorrCoeffFreqForecasts (y, z, zf, N, minNoObs, checkEnoughObs)

  % usage: [yhat, prms2, weights, R2, boundaries] = panelRandomCorrCoeffFreqForecasts (y, x, x_in_forecast_period, N, minNoObs, checkEnoughObs)
  %
  % forecasts from different forecasting methods for balanced panel with random
  % coefficients
  %
  % y - T*Nx1 vector of dependent variable
  % x - T*NxK, regressors without intercept
  % x_in_forecast_period - NxK, regressors in T+h without intercept
  % N - scaler, number of individuals
  % minNoObs - minimum number of observations for any individual, scalar
  % maxIterRE - max no iterations for BLUP RE estimation
  % checkEnoughObs - whether or not (1 or 0) to check if a unit has a
  %                  minimum no of observations
  %                  
  % yhat - Nx10, forecasts from different forecast methods
  % prms2 - parameter estimates
  % w - weight for combination forecast
  % R2 - cell with R2 for pooled and individual specific regressions
  
  % Andreas Pick 

  [NT,k] = size(z);
  T = NT/N;
  K = k+1;

  nMethods = 8;
  if nargin > 6
    if ~isempty(additionalMaterial{1})
      nMethods = nMethods + 2;
    end
    if length(additionalMaterial) > 1
      nMethods = nMethods + 2;
    end
  end

  yhat = nan(N,nMethods);

  % eliminating individuals without sufficient observations
  if checkEnoughObs == 1
    enoughObs = zeros(N,1);
    for ic = 1:N
      zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
      yi = y((ic-1)*T+1:ic*T);
      nansi = sum(isnan([zi yi]),2);
      zi = zi(nansi==0,:);
      if size(zi,1) < minNoObs
        continue
      elseif sum(isnan(zf(ic,:))) > 0
        continue
      elseif rank(zi'*zi) < K
        continue
      end
      enoughObs(ic) = 1;
    end
    z = z(kron(enoughObs,ones(T,1))==1,:);
    y = y(kron(enoughObs,ones(T,1))==1);
    zf = zf(enoughObs==1,:);
    N = sum(enoughObs);
    findObs = find(enoughObs==1);
  else
    findObs = 1:N;
  end

  % -----------------------------------------------------------------------
  % Forecasts:

  % Forecast #3: individual
  zz = zeros(K,K);
  zy = zeros(K,1);
  betaInd = nan(K,N);
  R2{2} = nan(N,1);
  SSE =  nan(N,1);
  SST = nan(N,1);
  Ti = nan(N,1);
  for ic = 1:N
    zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([zi yi]),2);
    Ti(ic) = sum(nansi==0);
    zi = zi(nansi==0,:);
    yi = yi(nansi==0);
    zz = zz + zi'*zi;
    zy = zy + zi'*yi;
    betaInd(:,ic) = pinv(zi'*zi)*(zi'*yi);

    % Forecast #3: Individual ==== 
    yhat(findObs(ic),3) = betaInd(1,ic) + zf(ic,:)*betaInd(2:end,ic);
    SSE(ic) = sum((yi - zi*betaInd(:,ic)).^2);
    SST(ic) = sum((yi - mean(yi)).^2);
    R2{2}(ic) = 1 - SSE(ic)/SST(ic);
  end

  % Forecast #1: Pooled ==== 
  betaPool = pinv(zz)*zy;
  yhat(findObs,1) = betaPool(1) + (zf*betaPool(2:end))';
  R2{1} = 1 - sum(SSE)/sum(SST);

  % Forecast #2: RE-BLUP ====
  w = [ones(NT,1) z];
  wf = [ones(N,1) zf];

  [yhat(findObs,2), ~] = GoldbergerREBLUP_FEinit (y, w, wf, N, T);

  % Forecast #4: FE ==== 
  [yhat(findObs,4), ~, betaFE] = FEforecastLargeN (y, z, zf, T);

  % Forecast #5: combination forecasts, individual and pooled =============
  
  [yhat_optW_pool, weight, boundPool, componentsPool] = optimalWeights_individual_pooled (yhat(findObs,3), yhat(findObs,1), y, w, wf, betaInd, betaPool, N);
  yhat(findObs,5) = yhat_optW_pool;

  % Forecast #6: combination forecasts, individual and FE
  
  [yhat_optW_FE, w_FE, boundFE, componentsFE] = optimalWeights_individual_FE (yhat(findObs,3), yhat(findObs,4), y, w, wf, betaInd, betaFE, N);
  yhat(findObs,6) = yhat_optW_FE;
  
  weights = [weight; w_FE];
  boundaries = [boundPool; boundFE];
  components = nan(2,2);
  components(1:2,1) = componentsPool;
  components(1:2,2) = componentsFE;

  % Forecast #7 and #8: equal weights
  weq = .5;
  yhat(findObs,7) = weq*yhat(findObs,3) + (1-weq)*yhat(findObs,1);
  yhat(findObs,8) = weq*yhat(findObs,3) + (1-weq)*yhat(findObs,4);

end
