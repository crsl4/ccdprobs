Integrate[
 Integrate[
  Integrate[
   x1^(a1 - 1)*x2^(a2 - 1)*
    x3^(a3 - 
       1)*(1 - x1 - x2 - x3)^(a4 - 1)/(lambda1*x1 + lambda2*x2 + 
        lambda3*x3 + lambda4*(1 - x1 - x2 - x3))^(a1 + a2 + a3 + 
        a4), {x1, 0, 1 - x2 - x3}, 
   Assumptions -> {Element[{x1, x2, x3, a1, a2, a3, a4, lambda1, 
       lambda2, lambda3, lambda4}, Reals], a1 > 0, a2 > 0, a3 > 0, 
     a4 > 0, lambda1 > 0, lambda2 > 0, lambda3 > 0, lambda4 > 0, 
     x1 > 0, x2 > 0, x3 > 0, x1 < 1, x2 < 1, x3 < 1}], {x2, 0, 
   1 - x3}], {x3, 0, 1}]