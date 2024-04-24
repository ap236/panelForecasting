function [yhat, w_FE, boundary, components] = optimalWeights_individual_FE (yhatInd, yhatFE, y, w, wf, betaInd, betaFEslope, N)

  % usage: [yhat, weight, boundary, components] = optimalWeights_individual_FE (yhatInd, yhatFE, y, w, wf, betaInd, betaFEslope, N)
  %
  % Calculates forecasts using optimal weights of Pesaran, Pick and
  % Timmermann (2024) that combine individual and FE forecasts.

  % Andreas Pick

  [NT, K] = size(w);
  k = K-1;
  T = NT/N;

  x = w(:,2:end);
  xf = wf(:,2:end);

  etaiFE = betaInd(2:end,:) - betaFEslope*ones(1,N);

  hNT = 0;   
  MT = eye(T) - ones(T,T)/T;

  DeltaFE = 0;
  
  for i = 1:N
    wi = w((i-1)*T+1:i*T,:);
    xi = x((i-1)*T+1:i*T,:);

    xfdmi = xf(i,:) - mean(xi);
    eihat = y((i-1)*T+1:i*T) - wi*betaInd(:,i); 
  
    HiT = xi'*MT*diag(eihat.^2)*MT*xi/T; 
    iQiT = inv(xi'*MT*xi/T);
    hNT = hNT + xfdmi*iQiT*HiT*iQiT*xfdmi'/N;

    DeltaFE = DeltaFE + xfdmi*etaiFE(:,i)*etaiFE(:,i)'*xfdmi'/N; 

  end
  
  w_FE = DeltaFE/(DeltaFE + hNT/T);

  boundary = 0;
  if w_FE < 0
    boundary = 1;
    w_FE = 0;
  elseif w_FE > 1
    boundary = 2;
    w_FE = 1;
  end

  components = [DeltaFE, hNT];

  yhat = w_FE*yhatInd + (1-w_FE)*yhatFE;

end