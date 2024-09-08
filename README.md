# ACS-MIP
Alternating Criteria Search is an heuristic algorithm that solves mixed-integer linear programming problems.

## Dependencies
- CPLEX 22.1.1

## Usage
Compile
```
gcc acs_mip.c -o acs_mip.out -lcplex -lm
```
Execute
```
Usage: ./acs_mip.out -i <intput> -o <output> -s <seed> -p <%varfixing>
   input: is a file with extension
      MPS, SAV, or LP (lower case is allowed)
   output: is a file with extension csv
   seed: any integer number
   %varfixing: is the percentage of the variable fixing
```
