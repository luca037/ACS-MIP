\* Example 2: Minimize cost with constraints and multiple integer variables *\

Minimize
 obj: 4 x1 + 3 x2 + 7 x3 + 5 x4
Subject To
 c1: 2 x1 + 3 x2 + x3 + 5 x4 <= 50
 c2: 4 x1 + x2 + 2 x3 + x4 <= 40
 c3: x1 + x2 + 3 x3 = 20
Bounds
 0 <= x1 <= 10
 0 <= x2 <= 15
 0 <= x3 <= 25
 5 <= x4 <= 12
Integers
 x1 x2 x3 x4
End
