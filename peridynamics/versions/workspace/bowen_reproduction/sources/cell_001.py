import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from dataclasses import dataclass, asdict
import json

np.set_printoptions(precision=5, suppress=True)
OUT = Path('bowen_stage1_outputs'); OUT.mkdir(exist_ok=True)
