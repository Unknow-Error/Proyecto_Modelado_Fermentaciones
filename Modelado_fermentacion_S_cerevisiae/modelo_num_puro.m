function [t_valores, conc_P_sol, conc_X_sol, conc_S_sol, conc_FN_sol, t_final_crecimiento] = modelo_num_puro (...
  conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_X_i, conc_S_i, conc_FN_i,...
  u_max, alpha, beta, n, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,...
  t_final_prueba, delta_t)
  % Función de funciones : Modela de forma diferencial la concentración del producto (etanol) en S.cerevisiae., su biomasa y sustrato.
  % Sistema de EDOs no lineales:
  % d[P]/dt = (alpha*umáx*(1-kp[P](t))^n + beta)*[x](t)
  % d[X]/dt = umax*(1-kp[P])^n*[X](t)
  % d[S]/dt = -(theta*umax*(1-kp[P])^n + lambda)*[x](t)
  % d[FN]/dt = -(c1/c3)*umax*(1-kp[P])^n *[x](t)
  % Parámetros cinéticos y del modelo : umáx, alpha, beta, n, conc_P_umbral, kp y conc_P_toxico.
  % Para la resolución numérica : parámetros t_final_prueba, delta_t, error_aceptable, max_iteracion, modo_exploratorio  y puntoEquilibrio véase "euler_implicito_ck.m".

  if nargin < 17, t_final_prueba = 50; endif
  if nargin < 18, delta_t = 0.1;  endif

  %Inicializar todas las salidas
  t_valores            = [];
  conc_P_sol           = [];
  conc_X_sol           = [];
  conc_S_sol           = [];
  conc_FN_sol          = [];
  t_final_crecimiento  = [];
  lambda = k_ms_2 + ms_1;
  theta = ((1/rend_v_y_x_s) + (k_ms_2/u_max));
  crecimiento_activo = true;

  % Las EDOs del sistema:
  d_conc_P_t = @(t, conc_P_t, conc_X_t) (alpha*monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n) + beta)*conc_X_t;
  d_conc_X_t = @(t, conc_P_t, conc_X_t) monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n)*conc_X_t;
  d_conc_S_t = @(t, conc_P_t, conc_X_t) -(theta*monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n) + lambda)*conc_X_t;
  d_conc_FN_t = @(t, conc_P_t, conc_X_t) -(c1/c3)*monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n)*conc_X_t;

  % Para ambos casos se puede modelar numéricamente hasta alcanzar 1/kp. Pueden ocurrir 2 cosas: O bien 1/kp = conc_P_toxico o bien 1/kp < conc_P_toxico.
  % En principio, kp afecta a la velocidad de crecimiento siendo 1/kp una concentración micoestática.
  % Mientras que conc_P_toxico se refiere a la concentración micocida de etanol (y que daría retroalimentación negativa final a la síntesis del mismo ya sea durante el crecimiento o en su actividad residual).
  % Otro detalle importante : Se formará etanol siempre y cuando exista sustrato. El sustrato se consume para crecimiento y para mantenimiento.
  % Si bien, está acoplado a la biomasa, si S -> 0 => Deja de crecer y deja de formarse producto, o bien, deja de formarse producto (si está ya en fase de mantenimiento).

  % Resolución mediante linealización local => Matriz Jacobiana:

  % Y = ([P], [X], [S], [FN]) => dY/dt = f(Y) = (f1([P] [X] [S] [FN]), f2([P] [X] [S] [FN]), f3([P] [X] [S] [FN]), f4([P] [X] [S] [FN]) = (d[P]/dt, d[X]/dt, d[S]/dt d[FN]/dt)
  % Supongamos un instante de tiempo

  % Valores de estado actuales (ejemplo)
  dY_dt = @(t, conc_P_t, conc_X_t)  [
          d_conc_P_t(t, conc_P_t, conc_X_t);
          d_conc_X_t(t, conc_P_t, conc_X_t);
          d_conc_S_t(t, conc_P_t, conc_X_t);
          d_conc_FN_t(t, conc_P_t, conc_X_t)
          ]; %Vector de Tasas.

  % Jij = dfi/dyj <- Matriz jacobiana
  % Para d[P]/dt : df1/d[P] = - (kp*alpha*umáx*n*(1-kp[P](t))^(n-1)*[x](t)
  %                df1/d[X] = alpha*umáx*(1-kp[P](t))^n + beta
  %                df1/d[S] = 0
  %                df1/d[FN] = 0
  % Para d[X]/dt : df2/d[P] = -(kp*umax*n*(1-kp[P])^(n-1))*[X](t)
  %                df2/d[X] = umax*(1-kp[P])^n
  %                df2/d[S] = 0
  %                df2/d[FN] = 0
  % Para d[S]/dt : df3/d[P] = +(kp*theta*umax*n*(1-kp[P])^(n-1))*[X](t)
  %                df3/d[X] = -(theta*umax*(1-kp[P])^n+lambda)
  %                df3/d[S] = 0
  %                df3/d[FN] = 0
  % Para d[FN]/dt : df4/d[P] = +(c1/c3)*kp*umax*n*(1-kp[P])^(n-1))*[X](t)
  %                 df4/d[X] = -(c1/c3)*umax*(1-kp[P])^(n)
  %                 df4/d[S] = 0
  %                 df4/d[FN] = 0

  %Valor inicial:
  Y_i            = [conc_P_i; conc_X_i; conc_S_i; conc_FN_i];
  valores_t      = 0:delta_t:t_final_prueba;
  N              = numel(valores_t);
  Y_aproximados  = zeros(4, N);
  dY_dt_aproximado = zeros(4, N);

  Y_aproximados(:,1)   = Y_i;
  dY_dt_aproximado(:,1)= dY_dt(valores_t(1), conc_P_i, conc_X_i);
  k_parada = N;

  for k = 1:N-1
    t_n   = valores_t(k);
    t_n_1  = valores_t(k+1);
    Y_n   = Y_aproximados(:,k);

    % Al agotarse FN, triggear bandera y truncar
    if (Y_n(1) > (1/kp) || Y_n(3) < 0 || Y_n(4) < 0)
      crecimiento_activo = false;
      Y_n(4) = 0;
    endif


    %Estimación de Y_n+1
    if ~crecimiento_activo
      % P evoluciona, X detenido, S solo mantenimiento, FN cero
      dY_dt_n = [
                beta*Y_n(2);
                0;
                -lambda * Y_n(2);
                0
      ];
    else
      dY_dt_n = dY_dt_aproximado(:,k);
    endif
    J_n = matriz_jacobiana(t_n, Y_n(1), Y_n(2), Y_n(3), Y_n(4), ...
          kp, conc_P_umbral, conc_P_toxico, u_max, alpha, beta, n, theta, lambda, c1, c3, crecimiento_activo);
    factor_jacobiano = eye(4)-delta_t*J_n;
    inversa_factor_jacobiano = inv(factor_jacobiano);
    Y_n_1 = Y_n + inversa_factor_jacobiano*delta_t*dY_dt_n;

    % Registrar t_final_crecimiento en el cruce de (1/kp) o si se agota el [S] o si se agota [FN]
    if (Y_n(1) < (1/kp) && Y_n_1(1) > (1/kp) && Y_n_1(3) > 0 && Y_n_1(4) > 0) || (Y_n_1(4) < 0 && Y_n_1(3) > 0) || (Y_n_1(3) < 0 && Y_n_1(4) > 0)
      t_final_crecimiento = t_n;
    endif

    % Condición de parada: S agotado o P tóxico
    if Y_n_1(3) <= 0 || Y_n_1(1) > conc_P_toxico
      k_parada = k;
      break;
    endif

    Y_aproximados(:,k+1)   = Y_n_1;

    % Estimar nuevo dY_dt
    J_n_1 = matriz_jacobiana(t_n_1, Y_n_1(1), Y_n_1(2), Y_n_1(3), Y_n_1(4), ...
           kp, conc_P_umbral, conc_P_toxico, u_max, alpha, beta, n, theta, lambda, c1, c3, crecimiento_activo);
    dY_dt_aproximado(:,k+1) = J_n_1*(Y_n_1-Y_n) + dY_dt_n;
  endfor

  t_valores = valores_t(1:k_parada)';
  conc_P_sol = Y_aproximados(1,1:k_parada)';
  conc_X_sol = Y_aproximados(2,1:k_parada)';
  conc_S_sol = Y_aproximados(3,1:k_parada)';
  conc_FN_sol = Y_aproximados(4, 1:k_parada)';
  if isempty(t_final_crecimiento)
    t_final_crecimiento = valores_t(k_parada);
  endif
endfunction

function matriz = matriz_jacobiana(...
  t, conc_P_t, conc_X_t, conc_S_t, conc_FN_t, kp, conc_P_umbral, conc_P_toxico,...
  u_max, alpha, beta, n, theta, lambda, c1, c3, crecimiento_activo)

  if (crecimiento_activo)
    matriz = [ -(kp)*alpha*u_max*n*max(0,((1-kp*conc_P_t)^(n-1)))*conc_X_t, alpha*max(0,monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n)) + beta, 0, 0;
               -(kp)*u_max*n*max(0,((1-kp*conc_P_t)^(n-1)))*conc_X_t, max(0,monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n)) , 0, 0;
                (kp)*theta*u_max*n*max(0,((1-kp*conc_P_t)^(n-1)))*conc_X_t, -(theta*max(0,monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n))  + lambda), 0, 0;
                (c1/c3)*kp*u_max*n*max(0,((1-kp*conc_P_t)^(n-1)))*conc_X_t, -(c1/c3)*max(0,monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n)) , 0, 0];
  else
    if (conc_P_t < conc_P_toxico && conc_S_t > 0)
      matriz = [0, beta, 0, 0;
              0, 0, 0, 0;
              0, lambda, 0, 0;
              0, 0, 0, 0];
    else
      matriz = zeros(4);
    endif
  endif
endfunction
