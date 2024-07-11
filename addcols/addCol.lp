\ENCODING=ISO-8859-1
\Problem name: mip_example.lp

Maximize
 obj1: x1 + 10 x2
Subject To
 c1: x1 + 7 x2 + 9 x3 <= 17.5
 c2: x1 + 10 x3 >= 1
Bounds
      x1 >= 0
      x2 >= 0
 0 <= x3 <= 10
Generals
 x1  x2  x3 
End
