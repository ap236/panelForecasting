function dm = PanelDieboldMarianoPPP_matlab (z, hvec, N, w, s)

  % usage: dm = PanelDieboldMarianoPPP_matlab (z, hvec, N, w, s)
  %
  % Panel Diebold Mariano test from Pesaran, Pick and Pranovich (2013, JoE, Appendix 2)
  %
  % Input: z - T*Nxh matrix of losses, for example difference in squared forecast errors
  %            stacked vectors of forecast losses per individual
  %            if not a balanced panel then still stacked forecast losses per
  %            individual and N (below) different
  %        hvec - vector of forecast horizons
  %        N - number of individuals - balanced panel: scalar
  %            unbalanced panel: vector with indicators which observations go
  %            with which individual
  %        w - weight vector to aggregate countries (optional, default 1/N)
  %        s - bandwidth for Bartlett kernel for h>1 (optional, default h-1)
  %
  % Output: dm - hx1 vector of panel DM statistics

  % Andreas Pick 

  if nargin == 3 % no weighting input then 1/N
    w = ones(N,1)/N;
  end
  if nargin < 5 % no bandwidth
    s = -9;
  end
  if min(hvec) < 1
    error('forecast horizon < 1 not allowed')
  end

  if size(z,2) == length(hvec)
    if length(N) == 1 % balanced panel
      if sum(sum(isnan(z))) == 0
        T = size(z,1)/N;
        balanced = 1;
      else
        balanced = 2;
      end
    else
      balanced = 0;
    end

  end

  % Now we start:

  dm = nan(length(hvec),1);
  hctr = 0;
  for h = hvec % loop over forecast horizons
    hctr=hctr+1;

    % balanced panel ---------
    if balanced == 1

      % version 1: max out on marix algebra - suitable for N not too large
      zbvec = nan(N,1);
      s2 = nan(N,1);
      for i = 1:N
        zbvec(i) = mean(z((i-1)*T+1:i*T,hctr));
        s2(i) = var(z((i-1)*T+1:i*T,hctr));
        if h > 1 % add correction for autocorrelation
          sadd = 0;
          if s < 0
            sh = h-1;
          else; sh = s;
          end
          for si = 1:sh
            sadd = sadd + (1-si/(1+sh))*(z((i-1)*T+1+si:i*T,hctr) - zbvec(i))'*(z((i-1)*T+1:i*T-si,hctr)-zbvec(i));
          end
          s2(i) = s2(i) + 2/T*sadd;
        end
      end
      zbar = sum(w.*zbvec);
      S = diag(s2);
      dm(hctr) = zbar./sqrt(w'*S*w/T);

      % unbalanced panels ---------
    else

      if balanced == 2
        T = size(z,1)/N;
        N = kron((1:N)',ones(T,1));
      end

      noInd = length(unique(N));
      zbvec = nan(noInd,1);
      s2 = nan(noInd,1);
      Ti = nan(noInd,1);
      noObs = [];
      for i = 1:noInd
        iN = unique(N);
        iN = iN(i);
        if iN == 0
          continue
        end
        zi = z(N==iN,hctr);
        zi = zi(~isnan(zi));
        if length(zi) < 3
          noObs = [noObs; i];
          continue
        end
        zbvec(i) = mean(zi);
        Ti(i) = sum(N==iN);
        s2(i) = var(zi);
        if h > 1 % add correction for autocorrelation
          sadd = 0;
          if s < 0
            sh = h-1;
          else
            sh = s;
          end
          for si = 1:sh
            sadd = sadd + (1-si/(1+sh))*(zi(si+1:end) - zbvec(i))'*(zi(1:end-si)-zbvec(i));
            endr
            s2(i) = s2(i) + 2/Ti(i)*sadd;
          end
        end
        if ~isempty(noObs)
          zbar(noObs) = [];
          s2(noObs) = [];
        end        
        zbar = sum(w.*zbvec);
        S = diag(s2);
        dm(hctr) = zbar./sqrt(w'*S*w*sum(1/Ti));
      end

    end
  end %endfunction
