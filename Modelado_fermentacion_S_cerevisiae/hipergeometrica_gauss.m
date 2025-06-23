function resultado = hipergeometrica_gauss(a,b,c,z, error_aceptable)
  % Aproximación de la función 2F1 hipergeometrica de Gauss por su Serie de Taylor.
  % 2F1(a,b;c;z) = sum_{k=0}^∞ [ ( (a)_k * (b)_k ) / ( (c)_k * k! ) ] * z^k
  % z : valor de entrada (es un numero complejo)
  % error_aceptable : diferencia mínima entre términos consecutivos
  % a, b y c modulan la serie hipergeométrica. Sus valores determinan si la serie converge a funciones conocidas como log(x), arctan(x), etc.
  % Calcula el valor por la serie de Taylor en el caso convergente si |z|<1. Para su extensión en casos |z|>1 se usa la extensión analítica por Transformación de Pfaff:
  % Transformación de Pfaff: 2F1(a,b;c;z) = (1-z)^(-a) 2F1(a,c-b;c;w) con w = z/(z-1) donde se tiene |z|>1 pero |w|<1.

  % Término inicial k = 0: (a)_0 = (b)_0 = (c)_0 = 1, y 0! = 1 => term = 1
  termino = 1;
  suma_actual = termino;
  k = 1;
  max_iter = 7500;

  if c == 0 || (c < 0 && floor(c) == c)
    error("El parámetro c no puede ser cero o un entero negativo para la función hipergeométrica 2F1.");
  endif

  if abs(z) <= 1
    % Caso base: calcula la serie directamente
    resultado = hipergeometrica_gauss_serie(a,b,c,z,error_aceptable);
  else
    % Transformación de Pfaff para |z| > 1
    w = z / (z - 1);
    factor = (1 - z)^(-a);
    resultado = factor * hipergeometrica_gauss_serie(a, c - b, c, w, error_aceptable);
  endif

endfunction

function res = hipergeometrica_gauss_serie(a,b,c,z,error_aceptable)
  % Serie de Taylor para 2F1(a,b;c;z) (usada internamente)

  termino = 1;
  suma_actual = termino;
  k = 1;
  max_iter = 7500;

  while k < max_iter
    termino = termino * (a + k - 1) * (b + k - 1) * z / ((c + k - 1) * k);
    suma_nueva = suma_actual + termino;

    if abs(suma_nueva - suma_actual) <= error_aceptable
      break;
    endif

    suma_actual = suma_nueva;
    k = k + 1;
  endwhile

  if k == max_iter
    warning(sprintf("hipergeometrica_gauss_serie: No convergió después de %d iteraciones para z = %f", max_iter, z));
  endif

  res = suma_actual;
endfunction
