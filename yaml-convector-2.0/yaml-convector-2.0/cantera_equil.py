import sys, json
import cantera as ct

if len(sys.argv) != 6:
    print("Usage: cantera_equil.py YAML_FILE PHASE_NAME COMPOSITION T P", file=sys.stderr)
    sys.exit(1)

yaml_path, phase_name, comp, T, P = sys.argv[1:6]
T = float(T)
P = float(P)

sol = ct.Solution(yaml_path, phase_name)
sol.TPX = T, P, comp
sol.equilibrate('TP')

for sp, x in zip(sol.species_names, sol.X):
    print(f"{sp} {x}")
