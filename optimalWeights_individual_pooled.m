function [yhat, weightPool, boundary, components] = optimalWeights_individual_pooled (yhatInd, yhatPooled, y, w, wf, betaInd, betaPool, N)

  % usage: [yhat, weight, boundary, components] = optimalWeights_individual_pooled (yhatInd, yhatPooled, y, w, wf, betaInd, betaPool, N)
  %
  % Calculates forecasts using optimal weights of Pesaran, Pick and
  % Timmermann (2024) that combine individual and pooled forecasts.

  % Andreas Pick 

  [NT, K] = size(w);
  T = NT/N;

  etaiPool = betaInd - betaPool*ones(1,N);

  hNT = 0;

  DeltaPool = 0;

  for i = 1:N
    wi = w((i-1)*T+1:i*T,:);
    yi = y((i-1)*T+1:i*T);
    eihat = yi - wi*betaInd(:,i);
    HiT = wi'*diag(eihat.^2)*wi/T;
    iQiT = inv(wi'*wi/T);
    hNT = hNT + wf(i,:)*iQiT*HiT*iQiT*wf(i,:)'/N;

    DeltaPool = DeltaPool + wf(i,:)*etaiPool(:,i)*etaiPool(:,i)'*wf(i,:)'/N;
  end

  weightPool = DeltaPool/(DeltaPool + hNT/T);
  weightPool_biasadjusted = (DeltaPool - hNT/T)/(DeltaPool);

  boundary = 0;
  if weightPool < 0
    weightPool = 0;
    boundary = 1;
  elseif weightPool > 1
    weightPool = 1;
    boundary = 2;
  end

  components = [DeltaPool, hNT];

  yhat = weightPool*yhatInd + (1-weightPool)*yhatPooled;

end