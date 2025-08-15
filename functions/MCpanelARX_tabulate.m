function MCpanelARX_tabulate (dataPath, doSuppl)

% prints tables of MC results

% Andreas Pick - 2025

if doSuppl == 1 % hierarchical Bayes
  MCrep = 1000;
else
  MCrep = 10000;
end

Tvec = [20, 50, 100]';
Nvec = [100, 1000];

rho0v = nan(3,1);
rho0v(1:3)  = [0.775, 0.688, 0.486]';
arhov = [0 0.5 1]';
g0v = ones(4,1)*0.1;
sig2gammav = [0, 0.1, 0.2]';
sig2alphav = [0.5, 0.5, 1]';
corxetav = [0, 1]';

% nLen = length(Nvec);
% tLen = length(Tvec);
% rLen = length(rho0v);

if doSuppl == 0
  disp('Table 1: Monte Carlo results')
  methodList = [1 2 4 8 5:7];
  nMethods = length(methodList);
elseif doSuppl == 1 % hierarchical Bayes
  disp('Table S.2: Monte Carlo results including hierarchical Bayes forecasts')
  methodList = [1 2 4 8 9:11 5 6];
  nMethods = length(methodList);
elseif doSuppl == 2 % Oracle and equal weights
  disp('Table S.4: Monte Carlo results for equally weighted and Oracle forecasts')
  methodList = 9:12;
  nMethods = length(methodList);
end

msfeC = cell(nMethods,1);

if doSuppl == 0

  for N = Nvec

    if N == 1000
      disp('Table S.3: Monte Carlo results, N=1000')
    end

    for MCpart = 1:2
      if MCpart == 1
        disp('Conditional on \kappa_i=0')
      else
        disp('Conditional on \kappa_i\pm 1')
      end

      for corrxb = corxetav'

        iT = 0;
        for T = Tvec'
          iT = iT+1;
          qr = 0;
          iA = 0;
          for rho0 = rho0v'
            iA = iA+1;
            qr = qr+1;
            sig2gamma = sig2gammav(iA);
            g0 = g0v(iA);
            sig2alpha = sig2alphav(iA);

            % load Monte Carlo results:
            fileStryerr = [dataPath 'results/MonteCarlo/PPTMCARX_Reps_'];
            infile1 = [fileStryerr ...
              num2str(MCrep) '_N_' num2str(N) ...
              '_T_' num2str(T) '_rho0_' num2str(rho0) '_sig2alpha_' ...
              num2str(sig2alpha) '_cor_beta_x_' ...
              num2str(corrxb)  '_MCpart_' num2str(MCpart)];

            if exist("infile1") == 0
              disp(['missing file with results: ' infile1])
            else
              load(infile1, '-mat')
            end

            mctr = 0;
            for iM =  methodList
              mctr = mctr+1;
              msfeC{mctr}(iA,iT) = msfei(iM);
            end

          end % rho
        end % T

        if corrxb == 0
          disp('\rho_{\gamma x}=0')
        else
          disp('\rho_{\gamma x}=0.5')
        end

        printString = '%1.1f & %1.1f &';
        for m = 1:(nMethods-1)
          printString = [printString ' %1.3f & %1.3f & %1.3f &&'];
        end
        printString = [printString ' %1.3f & %1.3f & %1.3f\\\\'];

        for m = 1:length(arhov)
          outData = [arhov(m) sig2alphav(m)];
          for jj = 1:nMethods
            outData = [outData msfeC{jj}(m,:)];
          end
          outstr = sprintf(printString, outData);
          disp(outstr)
        end

      end
    end % corrxb
  end % N

else

  for MCpart = 1:2

    if MCpart == 1
      disp('Conditional on \kappa_i=0')
    else
      disp('Conditional on \kappa_i\pm 1')
    end

    for N = Nvec

      for corrxb = corxetav'
        iT = 0;
        for T = Tvec'
          iT = iT+1;
          qr = 0;
          iA = 0;
          for rho0 = rho0v'
            iA = iA+1;
            qr = qr+1;
            sig2gamma = sig2gammav(iA);
            g0 = g0v(iA);
            sig2alpha = sig2alphav(iA);

            % load Monte Carlo results:
            fileStryerr = [dataPath 'results/MonteCarlo/PPTMCARX_Reps_'];
            infile1 = [fileStryerr ...
              num2str(MCrep) '_N_' num2str(N) ...
              '_T_' num2str(T) '_rho0_' num2str(rho0) '_sig2alpha_' ...
              num2str(sig2alpha) '_cor_beta_x_' ...
              num2str(corrxb)  '_MCpart_' num2str(MCpart)];
            if doSuppl == 2
              infile1 = [infile1 '_wSuppl'];
            end
            if doSuppl == 1
              infile1 = [infile1 '_wBayes'];
            end

            if exist("infile1") == 0
              disp(['missing file with results: ' infile1])
            else
              load(infile1, '-mat')
            end

            mctr = 0;
            for iM =  methodList
              mctr = mctr+1;
              msfeC{mctr}(iA,iT) = msfei(iM);
            end

          end % rho
        end % T

        if corrxb == 0
          disp(['N=' num2str(N) ', \rho_{\gamma x}=0'])
        else
          disp(['N=' num2str(N) ', \rho_{\gamma x}=0.5'])
        end

        printString = '%1.1f & %1.1f &';
        for m = 1:(nMethods-1)
          printString = [printString ' %1.3f & %1.3f & %1.3f &&'];
        end
        printString = [printString ' %1.3f & %1.3f & %1.3f\\\\'];

        for m = 1:length(arhov)
          outData = [arhov(m) sig2alphav(m)];
          for jj = 1:nMethods
            outData = [outData msfeC{jj}(m,:)];
          end
          outstr = sprintf(printString, outData);
          disp(outstr)
        end

      end % corrxb
    end % N
  end
end