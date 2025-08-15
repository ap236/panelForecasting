function [yhat, wis0, wis1] = panelRandomCoeffForecastsOldWeights (y, z, zf, N, minNoObs,checkEnoughObs)

  % usage: yhat = panelRandomCoeffForecastsOldWeights (y, x, x_in_forecast_period, N, minNoObs,checkEnoughObs)
  %
  % forecasts using weights in the working paper (Pesaran, Pick, Timmermann, 2022, CEPR WP)
  %
  % y - T*Nx1
  % x - T*NxK, regressors without intercept
  % x_in_forecast_period - NxK, regressors in T+h without intercept
  % N - scaler, number of individuals
  % minNoObs - min number of obs per individual
  % checkEnoughObs - if 1 check if each unit as enough observations
  %
  % yhat - Nx2, forecasts from different forecast methods
  % wis0,wis1 - scalars: how often bias adjusted weights are outside 0,1

  % Andreas Pick - 2025

  [NT,K] = size(z);
  T = NT/N;
  K = K+1;

  nMethods = 2;

  yhat = nan(N,nMethods);
  wis0 = 0;
  wis1 = 0;

  sigChoice = 2;

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

  % Forecasts:
  % -----------------------------------------------------------

  % 4) individual spec. OLS
  zz = zeros(K,K);
  zy = zeros(K,1);
  prms2 = nan(K,N);
  Ti = nan(N,1);
  yhatind = nan(N,1);
  for ic = 1:N
    zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([zi yi]),2);
    Ti(ic) = sum(nansi==0);
    zi = zi(nansi==0,:);
    yi = yi(nansi==0);
    zz = zz + zi'*zi;
    zy = zy + zi'*yi;
    prms2(:,ic) = pinv(zi'*zi)*(zi'*yi);
    %s2(ic) = 1/(Ti(ic)-K)*sum((yi-zi*prms2(:,ic)).^2);
    yhatind(ic) = prms2(1,ic) + zf(ic,:)*prms2(2:end,ic); %%% 3: individual
  end
  % pooled:
  prms1 = pinv(zz)*zy;
  yhatpool = prms1(1) + (zf*prms1(2:end))'; %%% 1: pool

  % 5) Mean group estimation forecast and rand coeff combination forecast
  wstaro = nan(N,1);
  wstaru = nan(N,1);
  wstnaive = nan(N,1);
  Omega_eta = cov(prms2');
  prmsMG = mean(prms2,2);
  sigmai2 = nan(N,1);
  sXX = zeros(K,K);
  for ic = 1:N

    zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([zi yi]),2);
    zi = zi(nansi==0,:);
    yi = yi(nansi==0);
    res = (yi - zi*prmsMG);
    if sigChoice == 1
      sigmai2(ic) = res'*res/(Ti(ic)+K); % mean group estimator of the variance
    elseif sigChoice == 2
      if ic == 1
        ExXXx = 0;
        for i = 1:N
          zi2 = [ones(T,1) z((i-1)*T+1:i*T,:)];
          yi = y((i-1)*T+1:i*T);
          nansi = sum(isnan([zi2 yi]),2);
          zi2 = zi2(nansi==0,:);
          if sum(isnan(zf(i,:))) > 0
            continue
          end
          ExXXx = ExXXx + 1/N*[1 zf(i,:)]*inv(zi2'*zi2)*[1 zf(i,:)]';
        end
      end
      sigmai2(ic) = res'*res/(Ti(ic)+ExXXx);
    end


    % 5) Naive combination weights

    ixixi = pinv(zi'*zi);
    xtp1 = [1 zf(ic,:)];
    wstaro(ic) = xtp1*Omega_eta*xtp1';
    wstaru(ic) = xtp1*(sigmai2(ic)*ixixi + Omega_eta)*xtp1';
    wstnaive(ic) = wstaro(ic)/wstaru(ic);
    yhat(findObs(ic),1) = (1-wstnaive(ic))*yhatpool(findObs(ic)) ...
      + wstnaive(ic)*yhatind(findObs(ic));

    % for first order optimal combination
    sXX = sXX + 1/N*sigmai2(ic)*pinv(zi'*zi);
  end



  % 10) Bias adjusted optimal weights

  Omega_etaC = Omega_eta - sXX; % expected value of Omega_eta
  for ic = 1:N
    zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([zi yi]),2);
    zi = zi(nansi==0,:);

    xtp1 = [1 zf(ic,:)];
    wstup = xtp1*(Omega_etaC)*xtp1';

    sXXi = sigmai2(ic)*pinv(zi'*zi);
    wstdown = xtp1*(Omega_etaC + sXXi)*xtp1';
    wstfo = wstup/wstdown;
    if wstfo > 1
      wstfo = 1;
      wis1 = wis1 + 1;
    elseif wstfo < 0
      wstfo = 0;
      wis0 = wis0 + 1;
    end
    yhat(findObs(ic),2) = (1-wstfo)*yhatpool(findObs(ic)) ...
      + wstfo*yhatind(findObs(ic));
  end

end
