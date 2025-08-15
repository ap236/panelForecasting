function [yhat, bRE] = GoldbergerREBLUP_FEinit (y, X, xf, N, T)

  % usage: yhat = GoldbergerREBLUP (y, X, xf, N, T)
  %
  % BLUP of random effects model

  % Andreas Pick - 2025

  if sum(X(:,1)==1)==size(X,1)
      X = X(:,2:end);
      xf = xf(:,2:end);
  end
  
  K = size(X,2);

  if N*T <= 10^6
    xdemean = X - kron(kron(eye(N),ones(1,T))*X/T,ones(T,1));
    ydemean = y - kron(kron(eye(N),ones(1,T))*y/T,ones(T,1));
  else
    xdemean = nan(N*T,K);
    ydemean = nan(N*T,1);
    for i = 1:N
      xdemean((i-1)*T+1:i*T,:) = X((i-1)*T+1:i*T,:) - ones(T,1)*mean(X((i-1)*T+1:i*T,:));
      ydemean((i-1)*T+1:i*T,:) = y((i-1)*T+1:i*T) - ones(T,1)*mean(y((i-1)*T+1:i*T));
    end
  end

  bFEslope = xdemean\ydemean;

  aFE = nan(N,1);
  sum1 = 0;
  sum2 = 0;
  MT = eye(T) - ones(T,T)/T;
  for i = 1:N
    yminxb = y((i-1)*T+1:i*T) - X((i-1)*T+1:i*T,:)*bFEslope;
    aFE(i) = mean(yminxb);
    sum1 = sum1 + yminxb'*MT*yminxb;
    sum2 = sum2 + (mean(y((i-1)*T+1:i*T)) - mean(X((i-1)*T+1:i*T,:))*bFEslope)^2;
  end

  sig2eps = 1/(N*(T-1)-K)*sum1;
  sig2alpha = 1/(N-K)*sum2 - sig2eps/T;
  if sig2alpha < 0
    sig2alpha = var(aFE);  
  end

  psi = sig2eps/(T*sig2alpha + sig2eps);

  iSig = 1/sig2eps*(MT + psi*(eye(T) - MT));
  XX = 0;
  Xy = 0;
  for i = 1:N
    XX = XX + [ones(T,1) X((i-1)*T+1:i*T,:)]'*iSig*[ones(T,1) X((i-1)*T+1:i*T,:)];
    Xy = Xy + [ones(T,1) X((i-1)*T+1:i*T,:)]'*iSig*y((i-1)*T+1:i*T);
  end
  bRE = inv(XX)*Xy;

  yhat1 = [ones(N,1) xf]*bRE; 
  uRE = y - [ones(N*T,1) X]*bRE;
  yhat21 = nan(N,1);
  for ii = 1:N
      yhat21(ii) = sum(uRE((ii-1)*T+1:ii*T));
  end
  yhat2 = sig2alpha/(T*sig2alpha + sig2eps)*yhat21;

  yhat = yhat1 + yhat2;

end
