function [x_valores, y_estimados] = euler_implicito_ck(x_inicial, x_final, y_inicial, ecuacion_diferencial, delta_x, error_aceptable, max_iteracion, modo_exploratorio, puntoEquilibrio)
  %Función que realiza Euler_implicito para una ecuación diferencial dada. Devuelve una matriz con valores de la función y para cada valor de x.
  % Método Euler implicito : Sea dy/dx = f(x,y) => (y_n+1-y_n)/h = f(x_n+1,y_n+1) => y_n+1 = h*f(x_n+1,y_n+1) + y_n. Tiene Error local de truncación O(h^2), error global N = (x1 - x0)/h : O(h)
  % Método Euler implicito con Crank-Nicolson (Trapecio implícito) : y_n+1 = y_n + (h/2)*(f(x_n,y_n) + f(x_n+1, y_n+1))
  % El anterior tiene un error de truncación O(h^3) , pero un error global de O(h^2).
  % La ecuación diferencial debe de ser de la forma lineal o no lineal del tipo @(x,y) o similar.
  % Por default tiene delta_x = 0.1, error_aceptable = 1e-8, max_iteracion = 5000 (para la convergencia), modo_exploratorio = false y puntoEquilibrio = [0,0]
  % Respecto al modo_exploratorio : si está en "true" permite desbloquear max_iteracion para saber cuantas iteraciones requieren realmente algunos valores para converger ante el error deseado.
  % Punto de equilibrio se refiere si hay un valor asintótico para la función determinando que se corte el computo una vez llegado allí dentro del error aceptable.
  % En puntoEquilibrio es un vector donde el primer valor corresponde a si se activa o no la comprobación y el segundo al valor de puntoEquilibrio final

  if nargin< 5, delta_x = 0.1;  endif
  if nargin< 6, error_aceptable = 1e-8; endif
  if nargin< 7, max_iteracion = 5000; endif
  if nargin< 8, modo_exploratorio= false; endif
  if nargin< 9, puntoEquilibrio= [0,0]; endif


  valores_x = x_inicial:delta_x:x_final;
  N = numel(valores_x);
  y_aproximado = zeros(1, N);
  y_aproximado(1) = y_inicial;

  salida_bucle = false;
  k_parada = N;

  for k = 1:N-1
    x_previo = valores_x(k);
    x_siguiente = valores_x(k+1);
    y_previo = y_aproximado(k);

    % Estimación inicial por Euler explícito para semilla.
    y_estimado = y_previo + delta_x * ecuacion_diferencial(valores_x(k), y_previo);

    % Iteración de punto fijo con criterio de error
    iteracion = 0;
    while iteracion < max_iteracion
      y_siguiente = y_previo + (delta_x/2) * (ecuacion_diferencial(x_previo, y_previo) + ecuacion_diferencial(x_siguiente, y_estimado));
      if abs(y_siguiente - y_estimado) <= error_aceptable
        y_estimado = y_siguiente;
        break
      endif
      y_estimado = y_siguiente;
      iteracion = iteracion + 1;
    endwhile

    %Chequeo de estabilización si el modo está activado
    if puntoEquilibrio(1) && y_estimado >= puntoEquilibrio(2)
      sprintf('euler_implicito_ck [puntoEquilibrio]: se terminó la ejecución al alcanzar el valor de equilibrio y=%.6f para x=%.4f', ...
              y_estimado, x_siguiente);
      salida_bucle = true;
      k_parada = k;
      break
    endif

    if iteracion == max_iteracion
      warning('euler_implicito_ck: no convergió en %d iteraciones para x=%.4f (y≈%.6f)', ...
              max_iteracion, x_siguiente, y_estimado);
      % Modo exploratorio : Intento casi ilimitado, con tope secundario
      if modo_exploratorio
        iteracion_exploratoria = 0;
        max_exploratorio = 5*max_iteracion;
        while iteracion_exploratoria < max_exploratorio
          y_siguiente = y_previo + (delta_x/2) * (ecuacion_diferencial(x_previo, y_previo) + ecuacion_diferencial(x_siguiente, y_estimado));
          if abs(y_siguiente - y_estimado) <= error_aceptable
            y_estimado = y_siguiente;
            break
          endif
          y_estimado = y_siguiente;
          iteracion_exploratoria = iteracion_exploratoria + 1;
        endwhile

        if iteracion_exploratoria == max_exploratorio
          warning('euler_implicito_ck [exploratorio]: No hubo convergencia tras %d totales para x=%.4f (y≈%.6f)', ...
                max_exploratorio, x_siguiente, y_estimado);
        else
          iteracion_total = iteracion_exploratoria + max_iteracion;
          warning('euler_implicito_ck [exploratorio]: convergió en %d iteraciones totales para x=%.4f (y≈%.6f)', ...
                iteracion_total, x_siguiente, y_estimado);
        endif
      endif
    endif

    y_aproximado(k+1) = y_estimado;
  endfor

  x_valores = valores_x(1:k_parada)';
  y_estimados = y_aproximado(1:k_parada)';

endfunction




