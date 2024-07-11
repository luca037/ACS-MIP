\ENCODING=ISO-8859-1
\Problem name: FMIP

Maximize
 obj1: dp1 + dp2 + dp3 + dn1 + dn2 + dn3
Subject To
 c1: x1 + 2 x2 + x3 + dp1 + dp2 + dp3 - dn1 - dn2 - dn3 <= 40
 c2: 2 x1 + x2 - x3 + dp1 + dp2 + dp3 - dn1 - dn2 - dn3 <= 30
 c3: x1 + x3 + dp1 + dp2 + dp3 - dn1 - dn2 - dn3  = 10
Bounds
      x1 = 20
 0 <= x2 <= 30
 0 <= x3 <= 10
Generals
 x1  x2  x3 
End
