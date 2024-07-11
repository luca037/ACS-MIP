\* Example 1: Maximize profit with constraints and integer variable *\

Maximize
 obj: 3 x1 + 5 x2 + 2 x3
Subject To
 c1: x1 + 2 x2 + x3 <= 40
 c2: 2 x1 + x2 - x3 <= 30
 c3: x1 + x3 = 10
Bounds
 0 <= x1 <= 20
 0 <= x2 <= 30
 0 <= x3 <= 10
Integers
 x2
End
