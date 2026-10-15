"""One-time auditable extraction; never executes notebook experiments."""
import ast,json
from pathlib import Path
root=Path(__file__).parent
nb=json.loads((root/'baselines/stage19/Bowen_RC_Chapter3_Stage19_Extended_Impact.ipynb').read_text())
out=root/'bowen_pd';out.mkdir(exist_ok=True)
(out/'__init__.py').write_text('')
groups={'geometry':[3,5,7,11,29,31,40,42,56,52], 'concrete':[94], 'rod':[44], 'interface':[74,75,96,118,104], 'contact':[54], 'legacy_runner':[120]}
header='import numpy as np\nfrom pathlib import Path\nfrom dataclasses import dataclass, asdict\nfrom numba import njit\nimport time, csv, json\n'
for name,cells in groups.items():
 parts=[header]
 for i in cells:
  tree=ast.parse(''.join(nb['cells'][i]['source']))
  if name=='geometry':
   # retain initialization; omit obsolete nonlinear concrete setup in cell 11
   if i==11: tree.body=[n for n in tree.body if isinstance(n,(ast.Assign,ast.FunctionDef))]
   # suppress only print calls
   tree.body=[n for n in tree.body if not (isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='print')]
  else: tree.body=[n for n in tree.body if isinstance(n,(ast.FunctionDef,ast.ClassDef))]
  parts+= [f'\n# Stage 19 notebook cell {i}; preserved mathematical operations\n',ast.unparse(tree)+'\n']
 (out/f'{name}.py').write_text(''.join(parts))
