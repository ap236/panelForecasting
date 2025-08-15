function getPanelR2Fun

% Andreas Pick - 2025

output = nan(3,4);

rho0v = [0.775, 0.688, 0.486];
g0 = [0.1 2/3 4/3]; % first number: mean of g0, second: multiple for first
                    % half, third: multiple for second half of units
sgv = [0,0.1,0.2];

	actr = 0;
	for abeta = [0,0.5,1]
		actr = actr+1;
		sg2 = sgv(actr);
		rho0 = rho0v(actr);
		
    output(actr,1:2) = [abeta rho0];
		for rbx = [0,0.5]
			
			output(actr,3+(rbx==.5)) = panelPR2 (rho0, abeta, g0, sg2, rbx);
			
		end
  end

disp('--- Table S.1 --- ')
printStr = '%1.3f %1.3f %1.3f %1.3f\n';
disp(sprintf(printStr,output'))