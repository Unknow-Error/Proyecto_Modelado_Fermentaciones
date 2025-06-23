function concentracion_FN = conc_FN_analitica_t(conc_P_t, conc_P_i, conc_FN_i, alpha, beta, c1, c3, u_max, kp, n)
  %Solución general analítica de la concentración de fuente de nidrógeno como [FN](t) = f([P](t)) del siguiente sistema de ecuaciones acoplado no lineal:
  %d_conc_P_t = alpha*d_conc_x_t + beta*conc_X_t;
  %d_conc_FN_t = -(c1/c3)*u_t*conc_X_t;

  if conc_P_t <= (1/kp)
    conc_FN_t = conc_FN_i + (c1/(c3*alpha*kp))*(kp*(conc_P_i-conc_P_t) ...
              -((1-kp*conc_P_t)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_t)^n, 10^(-8)) ...
              - (1-kp*conc_P_i)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_i)^n, 10^(-8))));
  else
    %Concentración [FN] mínima residual cuando ya creció toda la biomasa.
    conc_FN_t = conc_FN_i + (c1/(c3*alpha*kp))*((kp*conc_P_i-1)...
                 + (1-kp*conc_P_i)*hipergeometrica_gauss(1,1/n,(n+1)/n,-(alpha*u_max/beta)*(1-kp*conc_P_i)^n, 10^(-8)));
  endif

  concentracion_FN = conc_FN_t;

endfunction
