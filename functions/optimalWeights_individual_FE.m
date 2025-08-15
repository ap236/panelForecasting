function yhat = optimalWeights_individual_FE (yhatInd, yhatFE, y, w, wf, betaInd, betaFE, N) 

  % usage: [yhat, weight, boundary, components] = optimalWeights_individual_FE (yhatInd, yhatFE, y, w, wf, betaInd, alphaFE, betaFE, N)
  %
  % Calculates forecasts using optimal weights of Pesaran, Pick and
  % Timmermann (2025) that combine individual and FE forecasts.

  % Andreas Pick - 2025

  [NT, K] = size(w);
  k = K-1;
  T = NT/N;
  T = NT/N;
  if ~rem(T,2) % is even
    Thalf = T/2;
    Tsecondhalf = Thalf + 1;
  else
    Thalf = (T-1)/2;
    Tsecondhalf = Thalf + 2;
  end

  x = w(:,2:end);
  xf = wf(:,2:end);

  etaiMG = betaInd(2:end,:) - mean(betaInd(2:end,:),2)*ones(1,N);
  etaiFE = betaInd(2:end,:) - betaFE*ones(1,N);

  qNT = zeros(k,1);
  QNT = zeros(k,k);
  hNT = 0;   
  MT = eye(T) - ones(T,T)/T;
  psi1 = 0;
  psi2 = 0;
  DFE = 0;
  sigma2i = nan(N,1);
   
  for i = 1:N
    wi = w((i-1)*T+1:i*T,:);
    xi = x((i-1)*T+1:i*T,:);
    yi = y((i-1)*T+1:i*T);
    % half jackknife estimates
    ya = yi(1:Thalf);
    yb = yi(Tsecondhalf:end);
    wa = wi(1:Thalf,:);
    wb = wi(Tsecondhalf:end,:);
    ba = wa\ya;
    bb = wb\yb;
    bhjk = 2*betaInd(:,i) - 1/2*(ba+bb);
    
    xfdmi = xf(i,:) - mean(xi);
    eihat = yi - wi*betaInd(:,i); 
    sigma2i(i) = eihat'*eihat/(T-K);

    % estimates of h_NT
    HiT = sigma2i(i)*xi'*MT*xi/T; 
    iQiT = inv(xi'*MT*xi/T);
    hNT = hNT + xfdmi*iQiT*HiT*iQiT*xfdmi'/N;

    qNT = qNT + (xi'*MT*xi)*etaiMG(:,i)/(N*T);
    QNT = QNT + xi'*MT*xi/(N*T);

    % estimates of Delta
    DFE = DFE + xfdmi*etaiFE(:,i)*etaiFE(:,i)'*xfdmi'/N; 
    % estimate of psi
    psi1 = psi1 + T/N*(betaInd(2:end,i) - bhjk(2:end))'*(xfdmi'*xfdmi);
    psi2 = psi2 + T/N*(betaInd(2:end,i) - bhjk(2:end))'*(xfdmi'*xfdmi)*etaiMG(:,i);

  end
  
  psi = psi1*inv(QNT)*qNT - psi2;
  weight = (DFE - psi/T)/(DFE + hNT/T - 2*psi/T);

  if weight < 0
    weight = 0;
  elseif weight > 1
    weight = 1;
  end

  yhat = weight*yhatInd + (1-weight)*yhatFE;

end