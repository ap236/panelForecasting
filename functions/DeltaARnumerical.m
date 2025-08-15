function [Ey2eta, DeltaAR] = DeltaARnumerical (beta0)

	% numerically calculates values for Delta_AR for the 
	% appendix of Pesaran, Pick, Timmermann (2025)
	% for a = 1 and sigma^2 = 1

	% Andreas Pick - 2025

  rng(12345,'twister');
	etai = rand(10000,1)-.5;
	part1 = mean(etai.^2./(1-(beta0+etai).^2));
	part3 = mean(etai./(1-(beta0+etai).^2));
	part2 = 1/2*(log((1+beta0+1/2)/(1+beta0-1/2)) - log((1-beta0-1/2)/(1-beta0+1/2)));
	DeltaAR = (part1*part2 - part3^2)/part2;

  Ey2eta = 1/2*(-(1-beta0)*log((1-beta0-1/2)/(1-beta0+1/2)) ...
    -(1+beta0)*log((1+beta0+1/2)/(1+beta0-1/2)));

end