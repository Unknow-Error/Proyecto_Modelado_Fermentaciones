# Proyecto_Modelado_Fermentaciones

Este proyecto es un código implementado en MATLAB para la simulación de un proceso fermentativo anaeróbio de Saccharomyces cerevisiae mediante un sistema de ecuaciones diferenciales no lineales obtenidas a partir : 
1. Balance de materia en un sistema de cultivo Batch.
2. Modelo de Monod considerando inhibición por producto generalizada para la velocidad específica de crecimiento.
3. Modelo de Ludeking-Piret para velocidad de formación de etanol.
4. Modelo de Caja Negra para el balanceo de materia de la ecuación de crecimiento.
5. Modelo de Pirt para el análisis del sustrato destinado a mantenimiento celular.

## Uso

```bash
matlab -batch "run('/ruta/completa/Modelado_fermentacion_S_cerevisiae/Modelado_S_cerevisiae.m')"
```

## Rutinas Detalladas

1. **conc_FN_analitica_t.m**, **conc_S_analitica_t.m**, **conc_X_analitica_t.m** : Rutinas que definen las funciones analíticas que describen la evolución temporal en función de [P](t) de las concentraciones de [FN], [S] y [X].
2. **hipergeometrica_gauss.m** : Rutina que define la función hipergeométrica de Gauss.
3. **monod_modelo_inhb.m**: Rutina que describe el modelo de Monod con inhibición por producto.
4. **euler_implicito_ck.m**: Rutina que define una función de funciones para la resolución numérica de ecuaciones diferenciales mediante el método de Euler implicito usando la corrección de Crank-Nicolson (trapecio implicito) :  y_n+1 = y_n + (h/2)*(f(x_n,y_n) + f(x_n+1, y_n+1)), con un error de truncación O(h^3) , pero un error global de O(h^2). Por default tiene delta_x = 0.1, error_aceptable = 1e-8, max_iteracion = 5000 (para la convergencia), modo_exploratorio = false y puntoEquilibrio = [0,0].  Respecto al modo_exploratorio : si está en "true" permite desbloquear max_iteracion para saber cuantas iteraciones requieren realmente algunos valores para converger ante el error deseado. Punto de equilibrio se refiere si hay un valor asintótico para la función determinando que se corte el computo una vez llegado allí dentro del error aceptable.
5. **modelo_num_analitico.m** : Rutina que ejecuta la resolución analítica del modelo de ecuaciones diferenciales ordinarias no lineales que describen la evolución temporal de [S], [X], [YAN] y [X] para una fermentacion anaerobia de S. cerevisae en cervezas.
6. **modelo_num_puro.m** : Rutina que resuelve el sistema anteriormente mencionado de forma numérica realizando una aproximación lineal mediante el empleo de matrices jacobianas locales. 
