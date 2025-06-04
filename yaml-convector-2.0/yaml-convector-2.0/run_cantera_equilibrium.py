import sys
import cantera as ct

yaml_path = sys.argv[1]
phase_name = sys.argv[2]
composition = sys.argv[3]
T = float(sys.argv[4])
P = float(sys.argv[5])

sol = ct.Solution(yaml_path, phase_name)
sol.TPX = T, P, composition
sol.equilibrate('TP')

print(sol.report())
