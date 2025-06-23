%%% Modelo crecimiento S. cerevisiae anaeróbico en sistema Batch

% 1) Constantes del sistema:
%Cepa usada de referencia S. cerevisiae ATCC 4126
% Masa C-molar del etanol 23,034 g/cmol
% Formula empírica CH_1.613O_0.557N_0.158P_0.012S_0.003K_0.022Mg_0.003Ca_0.001 con una masa de 26,201 g/cmolX
u_max = 0.30; %Velocidad máxima de crecimiento en (h^-1), de 0.35 a 0.5 según medio y condiciones de pH.
rend_y_x_s = 0.25; %Rendimiento [cmolX/cmolS]
alpha = 3.5945; %Parámetro de Ludekin-Piret [cmolP/cmolX], referencia 2.805 g P/g X, 2.770 g P/g X, 3.105 g P/g X, a 3.515 g P/g X
beta = 0.00819; %Parámetro de Ludekin-Piret [cmolP/h*cmolX], referencia de 0.0011 a 0.0133 (g P/g X. h).
rend_v_y_x_s = 0.5445; %Parámetro de Pirt de Rendimiento verdadero [cmolX/cmolS], referencia de 0.5 gX/gS
ms_1 = 0.00482; %Parámetro de Pirt de mantenimiento no asociado a crecimiento [cmolS/cmolX*h], referencia de 0.0019 a 0.0086 g S/g X.h o bien 0.5 mmol/gX*h
k_ms_2 = 0.00551; %Parámetro de Pirt interno de mantenimiento (ms2) asociado a crecimiento [cmolS/cmolX*h], de 0.006 gS/gX*h (con umax de 0,3 h-1)
kp = 0.20566; %Parámetro de inhibición por producto en [cmolP/L] , [P] inhitoria = 112 g/L =  4,8624 cmolP/L => Kp = 1/[P]inhibitoria
n = 1; %Exponente y curva de la inhibición por producto. Con n = 1 inhibición lineal.
conc_P_umbral = 0.1*(1/kp); %Concentración donde el etanol empieza a realizar su efecto inhibitorio en el crecimiento en [cmolP/L]. Se puede aplicar desde el 10% como ejemplo.
conc_P_toxico = 5.7402; %Se encuentra que no hay más producción de etanol a los 115 g/L en la cepa elegida
c1 = 0.158; %Es el subíndice de la biomasa
c3 = 0.36101; %Es el subíndice de los YANs con CH_1,736O_0.27076N_0.36101

% 2) Condiciones iniciales del sistema:
% Supongamos para una cerveza de 12°Plato (12 g de extracto por cada 100 g de mosto) que en 1 L de agua sería 120 g/L
% El 80% (96 g/L) son azúcares fermentables siendo 8-12 g/L de Glucosa (10-15%), 0.8-1.6 g/L de Fructosa (1-2%), 40-48 g/L de maltosa (50-60%), 12-16 g/L de maltotriosa (15-20%) y <2g/L de sucrosa.
% Supongamos 10% GLucosa (C6H12O6), 1% de Fructosa(C6H12O6), 59% de maltosa (C12H22O11) y 10% de maltotriosa (C18H32O16): C_9.54H_17.5O_8.75 => CH_1.8344O_0.9172 (Masa C-molar : 28,5343 g/cmol)
% Nivel típico de YAN (mgN /L) de 150 a 300. 70 a 90% es FAN (Free aminoacid nitrogen, exceptuando prolina) y 10 a 30 % es NH3/NH4+ (por desaminaciones espontáneas de los tratamientos térmicos o actividad biológica)
% La formula molar promedio de los aminoacidosd libres es C_5.35H_7.85O_1.45N_1.45. La formula ponderada con los porcentajes (75% para FAN y 25% de NH3): C_2.77H_4.81O_0.75N_1, siendo en Cmol : CH_1,736O_0.27076N_0.36101
% La masa C-molar de esta composición de YAN es 23,146g/cmolYAN

conc_X_i = 0.028625; %Concentración inicial de biomasa en [cmolX/L], entre 0.5 a 1 g/L de levadura seca para 12°P a 15°P
conc_S_i = 3.3644; %Concentración inicial de sustrato en [cmolS/L], con 96g/L fermentables
conc_YAN_i = 0.5; %Concentración inicial de Yeast-assimilable nitrogen en [cmolYAN/L] para 300 mgYAN/L (0.013 cmol/L) típica de cerveza o 0.5cmolYAN/L de prueba.
conc_P_i = 0; %Concentración inicial del producto a t = ti
t_i = 0; %Tiempo inicial.

%% 3) Acoplamiento Biomasa-etanol: Sistema de ecuaciones diferenciales no lineales.
%d_conc_P_t = alpha*d_conc_x_t + beta*conc_X_t;
%u_t = monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n);
%d_conc_X_t = u_t*conc_X_t;

%% 3.1) Evolución temporal de la biomasa, sustrato e YAN (expresiones analíticas) y producción etanol (analítica-numérica)

% 3.1.1) Resolución numérica-analítica de la EDO de d[P]/dt anterior:
% EDO que describe el comportamiento del producto: Función anónima d[P]/dt = f(t, [P](t))
%d_conc_P_t = @(t, conc_P_t) (alpha*monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n) + beta)*conc_X_t(conc_P_t);
%Se hace desde modelo_conc_P_num_analitico.m
%conc_S_t = @(conc_P_t, t) conc_S_analitica_t(conc_P_t, conc_P_i, conc_S_i, conc_X_t((1/kp)), alpha, beta, rend_v_y_x_s, k_ms_2, ms_1, u_max, kp, n, t_final_crec, t);
% Lo mismo para la concentración de YAN porque ¿O se limita por azucar fermentables o se limita por YANs?

%% Resolución numérica de [P](t) : Requiere de la función [X](t) por su no linealidad y de [S](t) como limitante final (si se consume todo el sustrato).
% Contiene evaluación analítica de S[t] con los valores simulados de [P](t)

% Se analizó 4 formas de inhibición:
% Forma de inhibición lineal (n = 1) : u = umax(1-kp[P]) (Debería de dar un efecto logarítmico en la inhibición).
% Forma de inhibición cuadrática (n = 2) : u = umax(1-kp[P])^2 (Debería de dar un efecto arcotangente en la inhibición).
% Forma de inhibición radical (n = 1/2) : u = umax(1-kp[P])^(1/2)
% Forma de inhibición cúbica (n = 3): u = umax(1-kp[P])^3

t_computo_i = tic;
[t_sol_analitico_1_2, conc_P_sol_num_analitico_1_2, conc_X_sol_analitico_1_2, conc_S_sol_analitico_1_2, conc_FN_sol_analitico_1_2, t_final_crec_analitico_1_2] = modelo_num_analitico (...
    conc_X_i, conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_S_i, conc_YAN_i,...
    u_max, alpha, beta, (1/2), rend_v_y_x_s, k_ms_2, ms_1, c1, c3,...
    600, 0.01, 1e-6, 5000, true, [1,(1/kp)]);
t_computo_analitico_1_2 = toc(t_computo_i); %Medir el tiempo que tarda de computo para comparar algoritmos.

t_computo_i = tic;
[t_sol_analitico_1, conc_P_sol_num_analitico_1, conc_X_sol_analitico_1, conc_S_sol_analitico_1, conc_FN_sol_analitico_1, t_final_crec_analitico_1] = modelo_num_analitico (...
    conc_X_i, conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_S_i, conc_YAN_i, ...
    u_max, alpha, beta, n, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,...
    600, 0.01, 1e-6, 5000, true, [1,(1/kp)]);
t_computo_analitico_1 = toc(t_computo_i);

t_computo_i = tic;
[t_sol_analitico_2, conc_P_sol_num_analitico_2, conc_X_sol_analitico_2, conc_S_sol_analitico_2, conc_FN_sol_analitico_2, t_final_crec_analitico_2] = modelo_num_analitico (...
    conc_X_i, conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_S_i, conc_YAN_i,...
    u_max, alpha, beta, 2, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,...
    600, 0.01, 1e-6, 5000, true, [1,(1/kp)]);
t_computo_analitico_2 = toc(t_computo_i);

t_computo_i = tic;
[t_sol_analitico_3, conc_P_sol_num_analitico_3, conc_X_sol_analitico_3, conc_S_sol_analitico_3, conc_FN_sol_analitico_3, t_final_crec_analitico_3] = modelo_num_analitico (...
    conc_X_i, conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_S_i, conc_YAN_i,...
    u_max, alpha, beta, 3, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,...
    600, 0.01, 1e-6, 5000, true, [1,(1/kp)]);
t_computo_analitico_3 = toc(t_computo_i);

% 3.1.2) Resolución numérica pura usando jacobiano para linealizar localmente un sistema no lineal => modelo_conc_P_num_puro.m

t_computo_i = tic;
[t_sol_num_1_2, conc_P_sol_num_1_2, conc_X_sol_num_1_2, conc_S_sol_num_1_2, conc_FN_sol_num_1_2, t_final_crec_num_1_2] = modelo_num_puro (...
    conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_X_i, conc_S_i, conc_YAN_i, ...
    u_max, alpha, beta, (1/2), rend_v_y_x_s, k_ms_2, ms_1, c1, c3,  ...
    600, 0.01);
t_computo_numerico_1_2 = toc(t_computo_i);

t_computo_i = tic;
[t_sol_num_1, conc_P_sol_num_1, conc_X_sol_num_1, conc_S_sol_num_1, conc_FN_sol_num_1, t_final_crec_num_1] = modelo_num_puro (...
    conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_X_i, conc_S_i, conc_YAN_i, ...
    u_max, alpha, beta, n, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,  ...
    600, 0.01);
t_computo_numerico_1 = toc(t_computo_i);

t_computo_i = tic;
[t_sol_num_2, conc_P_sol_num_2, conc_X_sol_num_2, conc_S_sol_num_2, conc_FN_sol_num_2, t_final_crec_num_2] = modelo_num_puro (...
    conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_X_i, conc_S_i, conc_YAN_i, ...
    u_max, alpha, beta, 2, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,  ...
    600, 0.01);
t_computo_numerico_2 = toc(t_computo_i);

t_computo_i = tic;
[t_sol_num_3, conc_P_sol_num_3, conc_X_sol_num_3, conc_S_sol_num_3, conc_FN_sol_num_3, t_final_crec_num_3] = modelo_num_puro (...
    conc_P_i, conc_P_umbral, kp, conc_P_toxico, conc_X_i, conc_S_i, conc_YAN_i, ...
    u_max, alpha, beta, 3, rend_v_y_x_s, k_ms_2, ms_1, c1, c3,  ...
    600, 0.01);
t_computo_numerico_3 = toc(t_computo_i);

% 3.1.3) Valores finales y computo
sprintf('Los valores de tiempo de computo para la solución analítica fueron para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y n = 3: %.3f\n', ...
  t_computo_analitico_1_2, t_computo_analitico_1, t_computo_analitico_2, t_computo_analitico_3)

sprintf('Los valores de tiempo de computo para la solución numérica  fueron para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y n = 3: %.3f\n', ...
  t_computo_numerico_1_2, t_computo_numerico_1, t_computo_numerico_2, t_computo_numerico_3)

sprintf('Los valores finales de biomasa para la solución analítica  fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolX/L)\n',...
  conc_X_sol_analitico_1_2(end), conc_X_sol_analitico_1(end), conc_X_sol_analitico_2(end), conc_X_sol_analitico_3(end))

sprintf('Los valores finales de biomasa para la solución numérica fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolX/L)\n',...
  conc_X_sol_num_1_2(end), conc_X_sol_num_1(end), conc_X_sol_num_2(end), conc_X_sol_num_3(end))

sprintf('Los valores finales de etanol para la solución analítica  fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolP/L)\n',...
  conc_P_sol_num_analitico_1_2(end), conc_P_sol_num_analitico_1(end), conc_P_sol_num_analitico_2(end), conc_P_sol_num_analitico_3(end))

sprintf('Los valores finales de etanol para la solución numérica fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolP/L)\n',...
  conc_P_sol_num_1_2(end), conc_P_sol_num_1(end), conc_P_sol_num_2(end), conc_P_sol_num_3(end))

sprintf('Los valores finales de sustrato para la solución analítica  fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolS/L)\n',...
  conc_S_sol_analitico_1_2(end), conc_S_sol_analitico_1(end), conc_S_sol_analitico_2(end), conc_S_sol_analitico_3(end))

sprintf('Los valores finales de sustrato para la solución numérica fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolS/L)\n',...
  conc_S_sol_num_1_2(end), conc_S_sol_num_1(end), conc_S_sol_num_2(end), conc_S_sol_num_3(end))

sprintf('Los valores finales de YAN para la solución analítica  fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolYAN/L)\n',...
  conc_FN_sol_analitico_1_2(end), conc_FN_sol_analitico_1(end), conc_FN_sol_analitico_2(end), conc_FN_sol_analitico_3(end))

sprintf('Los valores finales de YAN para la solución numérica fueron  para para n = 1/2 : %.3f, n = 1: %.3f, n = 2 :  %.3f y  n = 3: %.3f  (cmolYAN/L)\n',...
  conc_FN_sol_num_1_2(end), conc_FN_sol_num_1(end), conc_FN_sol_num_2(end), conc_FN_sol_num_3(end))

% 4) Gráfica de las curvas
figure('Units', 'normalized', 'Position', [0.15 0.15 0.75 0.70]);

var_colors = {'r', 'g', 'b', 'm'};   % [P], [X], [S], [FN]
exponentes = {'n = 0.5', 'n = 1', 'n = 2', 'n = 3'};

datos_t      = {t_sol_analitico_1_2, t_sol_analitico_1, t_sol_analitico_2, t_sol_analitico_3};
datos_P      = {conc_P_sol_num_analitico_1_2, conc_P_sol_num_analitico_1, conc_P_sol_num_analitico_2, conc_P_sol_num_analitico_3};
datos_X      = {conc_X_sol_analitico_1_2,      conc_X_sol_analitico_1,      conc_X_sol_analitico_2,      conc_X_sol_analitico_3};
datos_S      = {conc_S_sol_analitico_1_2,      conc_S_sol_analitico_1,      conc_S_sol_analitico_2,      conc_S_sol_analitico_3};
datos_YAN     = {conc_FN_sol_analitico_1_2,     conc_FN_sol_analitico_1,     conc_FN_sol_analitico_2,     conc_FN_sol_analitico_3};

datos_t_num  = {t_sol_num_1_2,      t_sol_num_1,    t_sol_num_2,    t_sol_num_3};
datos_P_num  = {conc_P_sol_num_1_2, conc_P_sol_num_1, conc_P_sol_num_2, conc_P_sol_num_3};
datos_X_num  = {conc_X_sol_num_1_2, conc_X_sol_num_1, conc_X_sol_num_2, conc_X_sol_num_3};
datos_S_num  = {conc_S_sol_num_1_2, conc_S_sol_num_1, conc_S_sol_num_2, conc_S_sol_num_3};
datos_YAN_num = {conc_FN_sol_num_1_2,conc_FN_sol_num_1,conc_FN_sol_num_2,conc_FN_sol_num_3};

datos_tfin_crec     = [t_final_crec_analitico_1_2, t_final_crec_analitico_1, t_final_crec_analitico_2, t_final_crec_analitico_3];
datos_tfin_crec_num = [t_final_crec_num_1_2,       t_final_crec_num_1,       t_final_crec_num_2,       t_final_crec_num_3];

firstHandles = [];

for i = 1:4
    subplot(2,2,i);
    hold on;

    % Analítico
    h1 = plot(datos_t{i}, datos_P{i}, [var_colors{1} '-'], 'LineWidth',1.5);
    h2 = plot(datos_t{i}, datos_X{i}, [var_colors{2} '-'], 'LineWidth',1.5);
    h3 = plot(datos_t{i}, datos_S{i}, [var_colors{3} '-'], 'LineWidth',1.5);
    h4 = plot(datos_t{i}, datos_YAN{i},[var_colors{4} '-'], 'LineWidth',1.5);

    % Numérico (colores suaves)
    h5 = plot(datos_t_num{i}, datos_P_num{i}, 'o', 'MarkerEdgeColor',[1 0.8 0.6], 'MarkerFaceColor',[1 0.8 0.6]); %pastel naranja
    h6 = plot(datos_t_num{i}, datos_X_num{i}, 's', 'MarkerEdgeColor',[0.6 1 0.8], 'MarkerFaceColor',[0.6 1 0.8]); %verde pastel
    h7 = plot(datos_t_num{i}, datos_S_num{i}, 'd', 'MarkerEdgeColor',[0.6 0.8 1], 'MarkerFaceColor',[0.6 0.8 1]); %celeste pastel
    h8 = plot(datos_t_num{i}, datos_YAN_num{i},'^', 'MarkerEdgeColor',[0.8 0.6 1], 'MarkerFaceColor',[0.8 0.6 1]); %lila pastel

    % Capturamos sólo en la iteración 1
    if isempty(firstHandles)
        firstHandles = [h1,h5,h2,h6,h3,h7,h4,h8];
    end

    % Líneas verticales y textos
    y_max_analitico = max([datos_P{i}; datos_X{i}; datos_S{i}; datos_FN{i}]);
    y_max_num = max([datos_P_num{i}; datos_X_num{i}; datos_S_num{i}; datos_YAN_num{i}]);

    % t_crec analítico
    plot([datos_tfin_crec(i), datos_tfin_crec(i)], [0,y_max_analitico], 'b--','LineWidth',1.5);
    text(datos_tfin_crec(i)+0.2, 0.95*y_max_analitico, sprintf('t_{crec-analítico}=%.2f h',datos_tfin_crec(i)), 'FontSize',10,'Color','b');

    % t_crec numérico
    plot([datos_tfin_crec_num(i), datos_tfin_crec_num(i)], [0,y_max_num], 'k--','LineWidth',1.5);
    text(datos_tfin_crec_num(i)+0.2, 0.85*y_max_analitico, sprintf('t_{crec-numérico}=%.2f h',datos_tfin_crec_num(i)), 'FontSize',10,'Color','k');

    % t_sist analítico
    t_sist_analitico = datos_t{i}(end);
    plot([t_sist_analitico,t_sist_analitico],[0,y_max_analitico],'b:','LineWidth',1.5);
    text(t_sist_analitico+0.2,0.90*y_max_analitico,sprintf('t_{sist-analítico}=%.2f h',t_sist_analitico),'FontSize',10,'Color','b');

    % t_sist numérico
    t_sistema_numerico = datos_t_num{i}(end);
    plot([t_sistema_numerico,t_sistema_numerico],[0,y_max_num],'k:','LineWidth',1.5);
    text(t_sistema_numerico+0.2,0.80*y_max_num,sprintf('t_{sist-numérico}=%.2f h',t_sistema_numerico),'FontSize',10,'Color','k');

    title(['Inhibición: ', exponentes{i}], 'FontSize',14,'FontWeight','bold');
    xlabel('Tiempo (h)','FontSize',13);
    ylabel('Concentración (cmol/L)','FontSize',13);
    grid on; box on;
    hold off;
end

annotation('textbox',[0.10,0.95,0.9,0.10],...
    'String','Crecimiento de biomasa y producción de etanol – Diferentes valores de n - Analítico vs Numérico - Para 0.5 cmolYAN/L inicial','HorizontalAlignment','center',...
    'FontSize',18,'FontWeight','bold','EdgeColor','none');
legLabels = {...
    '[P] Etanol (Analítico)','[X] Biomasa (Analítico)',...
    '[S] Sustrato (Analítico)','[YAN] Nitrógeno (Analítico)',...
    '[P] Etanol (Numérico)','[X] Biomasa (Numérico)',...
    '[S] Sustrato (Numérico)','[YAN] Nitrógeno (Numérico)'};
handle_legend = legend (legLabels);
set(handle_legend, 'position', [ 0.15,   -0.12,  0.20,   0.30]);
set(handle_legend, 'fontsize', 9);
set(handle_legend, 'box', 'off');
set(handle_legend, 'orientation', 'horizontal');

% Determinación de la diferencia entre solución numérica y analítica.

N_errores = numel(datos_t);

MAE = struct('P',zeros(N_errores,1), 'X',zeros(N_errores,1), 'S',zeros(N_errores,1), 'YAN',zeros(N_errores,1)); % media de error absoluto (MAE)
RMSE = struct('P',zeros(N_errores,1), 'X',zeros(N_errores,1), 'S',zeros(N_errores,1), 'YAN',zeros(N_errores,1)); % raíz de error cuadrático medio (RMSE)
Einf = struct('P',zeros(N_errores,1), 'X',zeros(N_errores,1), 'S',zeros(N_errores,1), 'YAN',zeros(N_errores,1)); %error máximo

for i = 1:N_errores
  t_ana = datos_t{i};
  t_num = datos_t_num{i};

  % --- P ---
  P_ana_i = interp1(t_ana, datos_P{i}, t_num, 'spline');
  P_num_i = datos_P_num{i};
  error_P = P_num_i - P_ana_i;
  MAE.P(i)   = mean(abs(error_P));
  RMSE.P(i)  = sqrt(mean(error_P.^2));
  Einf.P(i)  = max(abs(error_P));

  % --- X ---
  X_ana_i = interp1(t_ana, datos_X{i}, t_num, 'spline');
  X_num_i = datos_X_num{i};
  error_X = X_num_i - X_ana_i;
  MAE.X(i)   = mean(abs(error_X));
  RMSE.X(i)  = sqrt(mean(error_X.^2));
  Einf.X(i)  = max(abs(error_X));

  % --- S ---
  S_ana_i = interp1(t_ana, datos_S{i}, t_num, 'spline');
  S_num_i = datos_S_num{i};
  error_S = S_num_i - S_ana_i;
  MAE.S(i)   = mean(abs(error_S));
  RMSE.S(i)  = sqrt(mean(error_S.^2));
  Einf.S(i)  = max(abs(error_S));

  % --- YAN ---
  YAN_ana_i = interp1(t_ana, datos_YAN{i}, t_num, 'spline');
  YAN_num_i = datos_YAN_num{i};
  error_YAN = YAN_num_i - YAN_ana_i;
  MAE.YAN(i)   = mean(abs(error_YAN));
  RMSE.YAN(i)  = sqrt(mean(error_YAN.^2));
  Einf.YAN(i)  = max(abs(error_YAN));

  % print resultados para la i-ésima
  fprintf('--- Experimento %d ---\n', i);
  fprintf(' P:   MAE=%.3e, RMSE=%.3e, E_inf=%.3e\n', MAE.P(i),  RMSE.P(i),  Einf.P(i));
  fprintf(' X:   MAE=%.3e, RMSE=%.3e, E_inf=%.3e\n', MAE.X(i),  RMSE.X(i),  Einf.X(i));
  fprintf(' S:   MAE=%.3e, RMSE=%.3e, E_inf=%.3e\n', MAE.S(i),  RMSE.S(i),  Einf.S(i));
  fprintf(' YAN:  MAE=%.3e, RMSE=%.3e, E_inf=%.3e\n\n', MAE.YAN(i), RMSE.YAN(i), Einf.YAN(i));
end
