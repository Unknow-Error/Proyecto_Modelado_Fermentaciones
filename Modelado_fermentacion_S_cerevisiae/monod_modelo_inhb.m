function resultado = monod_modelo_inhb(u_max, conc_P_t, conc_P_umbral, kp, n)
  %Función que retorna la ecuación de Monod para un cultivo en Batch según si hay o no inhibición por producto superando cierta concentración umbral.
  % u_max = velocidad máxima de crecimiento específico para el microoganismo en el medio seleccionado
  % conc_P_t = Función de la concentración de producto respecto del tiempo.
  % conc_P_umbral = concentración mínima en la cual el etanol comienza a ejercer su efecto micostático
  % kp = parámetro de la inhibición : Concentración máxima etanol producida durante crecimiento
  % n = parámetro de la inhibición : si es lineal, cuadrática, hiperbólica, etc.

  funcion_u = 1;

  if conc_P_t <= conc_P_umbral
    funcion_u = u_max; %Se comporta exponencial como en un Batch típico
  else
    if conc_P_t < (1/kp)
      funcion_u = u_max*(1-kp*conc_P_t)^n; %Se comporta de forma hipergeométrica habiendo un comportamiento pseudoexponencial para luego tener puntos de inflexión por la inhibición.
    else
      funcion_u = 0; %Se detiene el crecimiento.
    endif
  endif

  resultado = funcion_u;
endfunction
