function MCpanelARX_wrapper (runningNumber, MCpart, doSuppl, dataPath)

  % Wrapper function for MCpanelARX.m

  % Andreas Pick - 2025

  % === initial values of DGP ===
  arx = 1;

  Nvec = [100, 1000]'; 
  Tvec = [20, 50, 100]'; 
  rho0v = nan(3,1);
  rho0v(1:3)  = [0.775, 0.688, 0.486]';
  arhov = [0 0.5 1]';
  rhoctr = (1:3)';

  % ARX:
  if arx == 1
    g0 = 0.1;
    sig2gammav = [0, 0.1, 0.2]';
    sig2alphav = [0.5, 0.5, 1]';
    corxetav = [0, 1]';
    bigMatrix = [ ...
      kron(Nvec,ones(18,1)) ...
      kron(kron(ones(2,1),Tvec),ones(6,1)) ...
      kron(kron(ones(6,1),rhoctr),ones(2,1)) ...
      kron(ones(18,1),corxetav)
      ];
      if runningNumber > 36
        disp(['runningNumber too large ' num2str(runningNumber)])
        return
      end
    corrxb = bigMatrix(runningNumber,4);
    rho0c = bigMatrix(runningNumber,3);
    rho0 = rho0v(rho0c);
    arho = arhov(rho0c);
    sig2gamma = sig2gammav(rho0c);
    sig2alpha = sig2alphav(rho0c);
  % AR:
  elseif arx == 0 % pure AR panel
    g0 = 0;
    sig2gammav = [0, 0.1, 0.2]';
    sig2alphav = [0.5, 0.5, 1]';
    bigMatrix = [ ...
      kron(Nvec,ones(9,1)) ...
      kron(kron(ones(2,1),Tvec),ones(3,1)) ...
      kron(ones(6,1),rhoctr)
      ];
      if runningNumber > 18
        disp(['runningNumber too large ' num2str(runningNumber)])
        return
      end
    rho0c = bigMatrix(runningNumber,3);
    rho0 = rho0v(rho0c);
    arho = arhov(rho0c);
    sig2gamma = sig2gammav(rho0c);
    sig2alpha = sig2alphav(rho0c);
    corrxb = 0;
  end

  if sig2alpha < 0
    return;
  end

  N = bigMatrix(runningNumber,1);
  T = bigMatrix(runningNumber,2);

  format bank;
  disp(['runNo: ' num2str(runningNumber) ' N ' num2str(N) ' T ' num2str(T) ' kappa ' num2str(MCpart)])
  format short;

  if doSuppl == 1 % hierarchical Bayesian
    MCrep = 1000;
    MCpanelARX (MCrep, N, T, rho0, arho, corrxb, g0, sig2alpha, sig2gamma, MCpart, doSuppl, dataPath);
  else
    MCrep = 10000;
    MCpanelARX (MCrep, N, T, rho0, arho, corrxb, g0, sig2alpha, sig2gamma, MCpart, doSuppl, dataPath);
  end

end
