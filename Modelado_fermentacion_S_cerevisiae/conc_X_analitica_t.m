function concentracion_X = conc_X_analitica_t(conc_P_t, conc_P_i, conc_X_i, alpha, beta, u_max, kp, n)
  %Solución general analítica de la biomasa como [X](t) = f([P](t)) del siguiente sistema de ecuaciones acoplado no lineal:
  %d_conc_P_t = alpha*d_conc_x_t + beta*conc_X_t;
  %u_t = monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n);
  %d_conc_X_t = u_t*conc_X_t;

  if conc_P_t <= (1/kp)
    conc_X_t = conc_X_i + ((conc_P_t-conc_P_i)/alpha) ...
              + (1/(alpha*kp))*((1-kp*conc_P_t)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_t)^n, 10^(-8)) ...
              - (1-kp*conc_P_i)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_i)^n, 10^(-8)));
  else
    conc_X_t = conc_X_i + [((1/kp)-conc_P_i)/(alpha)]...
              -  ((1-kp*conc_P_i)/(alpha*kp))*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_i)^n, 10^(-8));
  endif

  concentracion_X = conc_X_t;

endfunction

