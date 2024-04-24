% making density plots for the relative msfe

% Andreas Pick

close('all')

savePlot = 0; % save plots
do_cdf = 0; % choose pdf (0) or cdf (1) version of empirical distributions

% -- Leave below as is ---
nModels = 3;
nMethods = 8;
rollwin = 60;
nWindows = length(rollwin);
modelVec = [1 4 5];

msfei = cell(nWindows,nModels);
noObs = cell(nWindows,nModels);
for i = 1:nWindows
  for j = 1:nModels
     msfei{i,j} = load([savePath 'msfei_' num2str(rollwin(i)) '_' num2str(modelVec(j)) '.txt']); % N x nMethods
     noObs{i,j} = load([savePath 'noObs_' num2str(rollwin(i)) '_' num2str(modelVec(j)) '.txt']);
  end
end

R2cell = cell(2,2);
for i = 1:nWindows
  for j = 1:nModels
    R2cell{i,j} = nan(size(msfei{i,j}));
    for l = 1:nMethods
      R2cell{i,j}(l,:) = msfei{i,j}(l,:)./msfei{i,j}(3,:);
    end
  end
end

colors = {'#77AC30' '#7E2F8E' '#D95319' '#0072BD' '#A2142F' '#4DBEEE' '#400EEE'}; % '#4D0232' '#477777'};
lines = {'-','-.',':','-','--','-.',':'}; %,'-','--'};

titletext{1,1} = 'AR';
titletext{2,1} = 'AR-PC';
titletext{3,1} = 'AR-X';

cntr = 0;
for i = 1:1 %:nWindows
  for j = 1:3 %:nModels
    cntr = cntr+1;
    figure(cntr)
    kctr = 0;
    for k = [1 2 4:8] 
      kctr = kctr+1;
      data = R2cell{i,j}(k,:);
      weights = 1./(noObs{1,1}(k,:)');
      x = 0.5:.005:1.5;
      y  = densityplotcalcs (data, x, 0.04); 
      if do_cdf == 1
        y = cumsum(y)./sum(y);
        ylim([0,1])
      end
      xlim([0.5,1.5])
      plot(x,y,'Color',colors{kctr},'LineStyle',lines{kctr},'LineWidth',2)
      hold on
    end
    xlabel('Ratio of MSFE')
    if do_cdf == 0
      legend('Pool', 'RE', 'FE', 'Emp.Bayes', 'Hier.Bayes', 'Comb (pool)', 'Comb (FE)')
    else
      legend('Pool', 'RE', 'FE', 'Emp.Bayes', 'Hier.Bayes', 'Comb (pool)', 'Comb (FE)','Location','northwest')
    end
    title([titletext{j,1}])
    if savePlot == 1
      if do_cdf == 1
        version = 'cdf';
      else
        version = 'pdf';
      end
      saveas(cntr,[savePath 'CPI_' titletext{j,1} '_' version '.pdf'],'pdf') 
      saveas(cntr,[savePath 'CPI_' titletext{j,1} '_' version '.fig'],'fig') 
    end
  end
end
hold off
