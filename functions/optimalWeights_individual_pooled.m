function yhat = optimalWeights_individual_pooled (yhatInd, yhatPooled, y, w, wf, betaInd, betaPool, N)

  % usage: yhat = optimalWeights_individual_pooled (yhatInd, yhatPooled, y, w, wf, betaInd, betaPool, N)
  %
  % Calculates forecasts using optimal weights of Pesaran, Pick and
  % Timmermann (2025) that combine individual and pooled forecasts.

  % Andreas Pick - 2025

  [NT, K] = size(w);
  T = NT/N;
  % for half-jackknife estimator:
  if ~rem(T,2) % is even
    Thalf = T/2;
    Tsecondhalf = Thalf + 1;
  else
    Thalf = (T-1)/2;
    Tsecondhalf = Thalf + 2;
  end

  etai = betaInd - mean(betaInd,2)*ones(1,N);
  etaP = betaInd - betaPool*ones(1,N);

  qNT = zeros(K,1);
  QNT = zeros(K,K);
  hNT = 0;
  psi1 = 0;
  psi2 = 0;
  Delta_pool = 0;
  sig2hat = nan(N,1);

  for i = 1:N
    wi = w((i-1)*T+1:i*T,:);
    yi = y((i-1)*T+1:i*T);
    % half jackknife estimates
    ya = yi(1:Thalf);
    yb = yi(Tsecondhalf:end);
    wa = wi(1:Thalf,:);
    wb = wi(Tsecondhalf:end,:);
    ba = wa\ya;
    bb = wb\yb;
    bhjk = 2*betaInd(:,i) - 1/2*(ba+bb);
    
    % (1) psi
    psi1 = psi1 + T/N*(betaInd(:,i) - bhjk)'*(wf(i,:)'*wf(i,:));
    psi2 = psi2 + T/N*(betaInd(:,i) - bhjk)'*(wf(i,:)'*wf(i,:))*etai(:,i);   

    qNT  = qNT + wi'*wi*etai(:,i)/(N*T);
    QNT  = QNT + wi'*wi/(N*T);
    
    % (2) hNT
    eihat  = yi - wi*betaInd(:,i);
    sig2hat(i) = sum(eihat.^2)/(T-K);
    HiT    = sig2hat(i)*(wi'*wi)/T;
    iQiT   = inv(wi'*wi/T);
 
    hNT    = hNT + wf(i,:)*iQiT*HiT*iQiT*wf(i,:)'/N;
    
    % (3) Delta 
    Delta_pool = Delta_pool + wf(i,:)*etaP(:,i)*etaP(:,i)'*wf(i,:)'/N; 

  end

  psi = psi1*inv(QNT)*qNT - psi2;
  weight = (Delta_pool - psi/T)/(Delta_pool + hNT/T - 2*psi/T); % weights using pooled eta

  if weight < 0
    weight = 0;
  elseif weight > 1
    weight = 1;
  end

  yhat = weight*yhatInd + (1-weight)*yhatPooled;
  
end