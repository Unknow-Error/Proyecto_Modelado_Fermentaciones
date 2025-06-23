function [t_valores, conc_P_sol, conc_X_sol, conc_S_sol, conc_FN_sol, t_final_crecimiento] = modelo_num_analitico(...
  conc_X_i, conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_S_i, conc_FN_i, ...
  u_max, alpha, beta, n, rend_v_y_x_s, k_ms_2, ms_1, c1, c3, ...
  t_final_prueba, delta_t, error_aceptable, max_iteracion, modo_exploratorio, puntoEquilibrio)
  % Función de funciones : Modela de forma numérica y analítica de la concentración del producto (etanol) en S.cerevisiae.
  % Detemina analíticamente a [X], [S] y [FN].
  % EDO: d[P]/dt = (alpha*umáx*(1-kp[P](t))^n + beta)*[x]([P](t))
  % Parámetros cinéticos y del modelo : umáx, alpha, beta, n, conc_P_umbral, kp y conc_P_toxico.
  % Para la resolución numérica : parámetros t_final_prueba, delta_t, error_aceptable, max_iteracion, modo_exploratorio  y puntoEquilibrio véase "euler_implicito_ck.m".

  if nargin < 14, t_final_prueba = 50; endif
  if nargin < 15, delta_t = 0.1;  endif
  if nargin < 16, error_aceptable = 1e-8; endif
  if nargin < 17, max_iteracion = 5000; endif
  if nargin < 18, modo_exploratorio= true; endif
  if nargin < 19, puntoEquilibrio= [1,(1/kp)]; endif

  % Para ambos casos se puede modelar numéricamente hasta alcanzar 1/kp. Pueden ocurrir 2 cosas: O bien 1/kp = conc_P_toxico o bien 1/kp < conc_P_toxico.
  % En principio, kp afecta a la velocidad de crecimiento siendo 1/kp una concentración micoestática.
  % Mientras que conc_P_toxico se refiere a la concentración micocida de etanol (y que daría retroalimentación negativa final a la síntesis del mismo ya sea durante el crecimiento o en su actividad residual).
  % Otro detalle importante : Se formará etanol siempre y cuando exista sustrato. El sustrato se consume para crecimiento y para mantenimiento.
  % Si bien, está acoplado a la biomasa, si S -> 0 => Deja de crecer y deja de formarse producto, o bien, deja de formarse producto (si está ya en fase de mantenimiento).

  %1) funciones auxiliares

  conc_X_t     = @(P) conc_X_analitica_t(P, conc_P_i, conc_X_i, alpha, beta, u_max, kp, n);
  d_conc_P_t   = @(t,P) (alpha*monod_modelo_inhb(u_max,P,conc_P_umbral,kp,n) + beta)*conc_X_t(P);
  conc_FN_t    = @(P)   conc_FN_analitica_t(P, conc_P_i, conc_FN_i, alpha, beta, c1, c3, u_max, kp, n);
  conc_S_t = @(P, conc_X_final, t_final_crecimiento, t) conc_S_analitica_t(P, conc_P_i, conc_S_i, conc_X_final, ...
                alpha, beta, rend_v_y_x_s, k_ms_2, ms_1, u_max, kp, n, t_final_crecimiento, t);

  %2) fase numérica inicial
  [t_valores, conc_P_sol] = euler_implicito_ck(0, t_final_prueba, conc_P_i, d_conc_P_t, delta_t, error_aceptable, max_iteracion, modo_exploratorio, puntoEquilibrio);
  t_valores = t_valores(:);
  conc_P_sol = conc_P_sol(:);
  t_final_crecimiento = t_valores(end);
  conc_P_ultm_crec = conc_P_sol(end);; %En el caso de que no se limite ni por FN o S , se limita por 1/kp en principio o por el valor alcanzado por equilibrio en la resolucion anterior.
  conc_X_final = conc_X_t(conc_P_ultm_crec); %Para testeo


  %3) fase de crecimiento => Se detiene [X] si se agota [FN] o [S]
  [t_parciales, P_parcial, X_parcial, S_parcial, FN_parcial, t_corte] = ...
    truncar_por_S_o_FN(t_valores, conc_P_sol, conc_X_t, conc_X_final, conc_S_t, conc_FN_t, t_final_crecimiento);

  t_valores = t_parciales(:);
  conc_S_sol = S_parcial(:);
  conc_FN_sol= FN_parcial(:);
  conc_X_sol = X_parcial(:);
  conc_P_sol = P_parcial(:);
  conc_P_ultm_crec = conc_P_sol(end);
  if ~isempty(t_corte)
    t_final_crecimiento = t_corte; %Nuevo tiempo final de crecimiento si se agota S o FN
  end

  %4) fase post-crecimiento (aplica si hay S residual)
  if (conc_P_ultm_crec < conc_P_toxico && conc_S_sol(end) > 0)
    conc_X_final = conc_X_sol(end);
    t_final_sistema = (conc_P_toxico - conc_P_ultm_crec)/(beta*conc_X_final) + t_valores(end);
    tspan = (t_valores(end) + delta_t) : delta_t : t_final_sistema;
    if numel(tspan)>=2, tspan = tspan(2:end); endif
    tspan = tspan(:);

    %Calculo del sustrato y ver cuando termina
    S_post = conc_S_sol(end) - (k_ms_2+ms_1)*conc_X_final*(tspan - tspan(1));
    indice_final = find(S_post <= 0, 1, 'first');
    tspan = tspan(1:indice_final);
    S_post = S_post(1:indice_final);

    %Calcular el resto
    P_post = conc_P_ultm_crec + beta*conc_X_final*(tspan - tspan(1));
    FN_post = conc_FN_sol(end) * ones(indice_final,1);
    X_post = conc_X_final* ones(indice_final,1);

    t_valores   = [t_valores; tspan];
    conc_P_sol  = [conc_P_sol; P_post(:)];
    conc_S_sol  = [conc_S_sol; S_post(:)];
    conc_FN_sol = [conc_FN_sol; FN_post(:)];
    conc_X_sol  = [conc_X_sol; X_post(:)];
  endif

endfunction

% Función helper interno para truncar P y S cuando S llega a cero
function [t_salida, P_salida, X_salida, S_salida, FN_salida, t_corte] = truncar_por_S_o_FN(t_entrada, P_entrada, conc_X_t, X_final, conc_S_t, conc_FN_t, t_final)
  if isempty(t_entrada)
    t_salida = [];
    P_salida = [];
    X_salida = [];
    S_salida = [];
    FN_salida= [];
    t_corte  = [];
    return
  endif

  N = numel(t_entrada);
  FN_vec = zeros(N,1);
  S_vec  = zeros(N,1);

  % Calcular FN y S, y verificar si se anula o se hace negativo:
  for i = 1:N
    P_val = P_entrada(i);
    t_val = t_entrada(i);

    % --- FN ---
    FN_i = conc_FN_t(P_val);
    if FN_i < 0
        FN_i = 0;
    end
    FN_vec(i) = FN_i;

    % --- S ---
    S_i = conc_S_t(P_val, X_final, t_final, t_val);
    if S_i < 0
        S_i = 0;
    end
    S_vec(i) = S_i;

    if S_i ==0 || FN_i == 0
      break;
    endif
  endfor

  indice_FN = find(FN_vec <= 0, 1, 'first');
  indice_S  = find(S_vec  <= 0, 1, 'first');

  % En el caso donde no se agoto ninguno todavía:
  if isempty(indice_S),  indice_S  = numel(t_entrada);  endif
  if isempty(indice_FN), indice_FN = numel(t_entrada);  endif

  % El crecimiento “acaba” en el primer evento (S o FN)
  indice_corte = min(indice_S, indice_FN);
  t_corte   = t_entrada(indice_corte);

  % Redimensionar los vectores al nuevo valor final de crecimiento => Salida o return de esta función durante la fase de crecimiento
  if ~isempty(indice_corte) && indice_corte > 1 && (t_entrada(indice_corte-1) < t_final)
    t_salida = t_entrada(1:indice_corte-1);
    P_salida = P_entrada(1:indice_corte-1);
    X_salida = arrayfun(@(P_val) conc_X_t(P_val), P_salida);
    S_salida = S_vec(1:indice_corte-1);
    FN_salida = FN_vec(1:indice_corte-1);
  endif
endfunction
