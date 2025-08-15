function R2 = panelPR2 (rho0, arho, g0, sg2, rbx)

  % Andreas Pick, 2025

  Vg = sg2*(rbx^2/(1-rbx^2)+1);
  if isscalar(g0)
    Eg2 = g0^2 + Vg;
  else
    g02 = g0(1)^2 * (g0(2)^2/2 + g0(3)^2/2);
    Eg2 = g02 + Vg;
  end
  if arho == 0
     Er2 = 1/(1-rho0^2);
  else
     Er2 = 1/(2*arho)*(log((1+rho0+arho/2)/(1+rho0-arho/2)) - ...
        log((1-rho0-arho/2)/(1-rho0+arho/2)));
  end
  R2 = (Eg2*Er2 + Er2 - 1)/(Eg2*Er2 + Er2);

end
