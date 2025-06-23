function concentracion_S = conc_S_analitica_t(conc_P_t, conc_P_i, conc_S_i, conc_X_max, alpha, beta, y_xs_verd, k_ms2, ms1, u_max, kp, n, tbiof, t)
  %Solución general analítica de la biomasa como [X](t) = f([P](t)) del siguiente sistema de ecuaciones acoplado no lineal:
  %d_conc_P_t = alpha*d_conc_x_t + beta*conc_X_t;
  %Pirt: 1/yxs = 1/yxs_v + ms/umax, con ms = ms1 + ms2 = ms1 + k(1-u_t/u_max)
  %u_t = monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n);
  %d_conc_S_t = -(theta*u_t+lambda)*conc_X_t;

  lambda = k_ms2 + ms1;
  theta = ((1/y_xs_verd) + (k_ms2/u_max));

  if conc_P_t <= (1/kp)
    conc_S_t = conc_S_i + (theta/alpha)*(conc_P_i - conc_P_t) ...
              + ((alpha*lambda-beta*theta)/(alpha*beta*kp))*(((1-kp*conc_P_t)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_t)^n, 10^(-8)) ...
              - (1-kp*conc_P_i)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_i)^n, 10^(-8))));
  else
    %Concentración [S] cuando ya creció toda la biomasa.
    conc_S_tbiof = conc_S_i + (theta/alpha)*(conc_P_i-(1/kp))...
                 - ((alpha*lambda-beta*theta)/(alpha*beta*kp))*((1-kp*conc_P_i)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_i)^n, 10^(-8)));
    %Consumo por mantenimiento.
    conc_S_t = conc_S_tbiof - lambda*conc_X_max*(t-tbiof);
  endif

  concentracion_S = conc_S_t;

endfunction
