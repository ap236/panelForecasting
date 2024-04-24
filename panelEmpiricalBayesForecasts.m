function yhat = panelEmpiricalBayesForecasts (y, z, zf, N, minNoObs, checkEnoughObs)

  % usage: yhat = panelEmpiricalBayesForecasts (y, z, zf, N, minNoObs, checkEnoughObs)
  %
  % forecasts empirical Bayes (Hsiao et al 1999)
  %
  % y - T*Nx1 vector of dependent variable
  % z - T*NxK, regressors without intercept
  % zf - NxK, regressors in T+h without intercept
  % N - scaler, number of individuals
  % minNoObs - minimum number of observations for any individual, scalar
  % checkEnoughObs - whether or not (1 or 0) to check if a unit has a
  %                  minimum no of observations
  %
  % yhat - Nx1, forecasts from different forecast methods

  % Andreas Pick 

  [NT,K] = size(z);
  T = NT/N;
  K = K+1;

  yhat = nan(N,1);

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

  % version with unadjusted variance
  betai = nan(K,N);
  sig2i = nan(N,1);
  for ic = 1:N
    zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([zi yi]),2);
    zi = zi(nansi==0,:);
    yi = yi(nansi==0);
    betai(:,ic) = inv(zi'*zi)*(zi'*yi);
    eps = yi - zi*betai(:,ic);
    sig2i(ic) = eps'*eps/(size(eps,1)-K);
  end

  barbeta = mean(betai,2);

  Omega = zeros(K,K);
  for i = 1:N
    Omega = Omega + (betai(:,i)-barbeta)*(betai(:,i)-barbeta)'/N;
  end
  iOmega = inv(Omega);

  betaEB = nan(K,N);
  for ic = 1:N
    zi = [ones(T,1) z((ic-1)*T+1:ic*T,:)];
    yi = y((ic-1)*T+1:ic*T);
    nansi = sum(isnan([zi yi]),2);
    zi = zi(nansi==0,:);
    yi = yi(nansi==0);
    betaEB(:,ic) = inv(1/sig2i(ic)*zi'*zi + iOmega)*(zi'*yi/sig2i(ic) + iOmega*barbeta);
    yhat(findObs(ic)) = betaEB(1,ic) + zf(ic,:)*betaEB(2:end,ic);
  end


end
