# Hierarchical linear solver overhead for an OpenAeroStruct multipoint problem
This repo hosts a sample script that demonstrates the computational overhead where the OAS-level linear solver gets a lot of unnecessary calls with RHS vector = 0.
This happens when the outputs of the OAS analyses (e.g., CL) are concatenated into vector variables.pip

## Problem setup
- Problem: multi-point OpenAeroStruct
- Goal: compute the derivatives of `CL` of each point in the reverse mode.
- Linear solver: `PETScKrylov` for each aerostructural point.

## Package versions
- OpenMDAO 3.28.1.dev0
- OpenAeroStruct 2.7.0

## Runscript and results

```
python run_multipoint_RHS0_overhead.py
```
runs the 2-point OAS analyses and compute the total derivatives.

### Case 1: Constraining CL from each point separately
```python
# CL constraints
for i in range(num_points):
    prob.model.add_constraint('AS_point_' + str(i) + '.CL', equals=0.5, ref=0.1)
```
The code outputs are shown below.
There are 2 OAS adjoint solves; no overhead.
```
--- constraining CL of each point ---
LN: PETScKrylov 0 ; 0.000760291935 1
LN: PETScKrylov 1 ; 1.03147509e-07 0.000135668293
LN: PETScKrylov 2 ; 1.38659576e-09 1.82376755e-06
LN: PETScKrylov 3 ; 3.73692855e-12 4.91512323e-09
LN: PETScKrylov 0 ; 0.00164369345 1
LN: PETScKrylov 1 ; 1.00038291e-07 6.08618907e-05
LN: PETScKrylov 2 ; 4.69893243e-09 2.85876447e-06
LN: PETScKrylov 3 ; 1.01143174e-11 6.15340858e-09
```

### Case 2: Constraining CL vector
```python
# CL constraints
prob.model.add_constraint('CL_vec', equals=0.5, ref=0.1)  # CL_vec has a length of num_points
```
The concatenations of `CL` is done by the `Scalars2Vector` component implemented in `run_multipoint_RHS0_overhead.py`.  
As a result, there are 4 OAS adjoint solvers: 2 with RHS=0 and 2 with non-zero RHS.
```
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.000760291935 1
LN: PETScKrylov 1 ; 1.03147509e-07 0.000135668293
LN: PETScKrylov 2 ; 1.38659576e-09 1.82376755e-06
LN: PETScKrylov 3 ; 3.73692855e-12 4.91512323e-09
LN: PETScKrylov 0 ; 0.00164369345 1
LN: PETScKrylov 1 ; 1.00038291e-07 6.08618907e-05
LN: PETScKrylov 2 ; 4.69893243e-09 2.85876447e-06
LN: PETScKrylov 3 ; 1.01143174e-11 6.15340858e-09
LN: PETScKrylov 0 ; 0 0
```

### Case 2 with 10 analysis points
There are `(n**2 - n)` RHS=0 calls, where `n=10` for 10-point problem.
The overhead becomes significant for large `n`.

```
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.000760291935 1
LN: PETScKrylov 1 ; 1.03147509e-07 0.000135668293
LN: PETScKrylov 2 ; 1.38659576e-09 1.82376755e-06
LN: PETScKrylov 3 ; 3.73692855e-12 4.91512323e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.000858888928 1
LN: PETScKrylov 1 ; 1.02704266e-07 0.000119578053
LN: PETScKrylov 2 ; 1.7729998e-09 2.06429463e-06
LN: PETScKrylov 3 ; 3.90145543e-12 4.54244467e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.000957396489 1
LN: PETScKrylov 1 ; 1.02282846e-07 0.000106834365
LN: PETScKrylov 2 ; 2.15786071e-09 2.25388409e-06
LN: PETScKrylov 3 ; 4.16797038e-12 4.35344231e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.00105580562 1
LN: PETScKrylov 1 ; 1.01884303e-07 9.64991103e-05
LN: PETScKrylov 2 ; 2.53967301e-09 2.40543615e-06
LN: PETScKrylov 3 ; 4.56046071e-12 4.31941316e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.00115410736 1
LN: PETScKrylov 1 ; 1.01509741e-07 8.79551975e-05
LN: PETScKrylov 2 ; 2.91714937e-09 2.52762391e-06
LN: PETScKrylov 3 ; 5.09630769e-12 4.41580035e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.00125229277 1
LN: PETScKrylov 1 ; 1.0116031e-07 8.07800798e-05
LN: PETScKrylov 2 ; 3.28909184e-09 2.62645598e-06
LN: PETScKrylov 3 ; 5.78542042e-12 4.61986251e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.00135035294 1
LN: PETScKrylov 1 ; 1.00837202e-07 7.46747015e-05
LN: PETScKrylov 2 ; 3.65434815e-09 2.70621706e-06
LN: PETScKrylov 3 ; 6.63175436e-12 4.9111267e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.00144827899 1
LN: PETScKrylov 1 ; 1.0054165e-07 6.94214653e-05
LN: PETScKrylov 2 ; 4.01179717e-09 2.77004444e-06
LN: PETScKrylov 3 ; 7.63583257e-12 5.27234918e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.0015460621 1
LN: PETScKrylov 1 ; 1.00274916e-07 6.48582718e-05
LN: PETScKrylov 2 ; 4.36034599e-09 2.82029163e-06
LN: PETScKrylov 3 ; 8.79691118e-12 5.68988217e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0.00164369345 1
LN: PETScKrylov 1 ; 1.00038291e-07 6.08618907e-05
LN: PETScKrylov 2 ; 4.69893243e-09 2.85876447e-06
LN: PETScKrylov 3 ; 1.01143174e-11 6.15340858e-09
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
LN: PETScKrylov 0 ; 0 0
```



