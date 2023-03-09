import sys
from pathlib import Path

project_dir = Path(__file__).parents[2]
cascade_dir = project_dir / "cascade"
sys.path.append(str(project_dir))
sys.path.append(str(cascade_dir))
# # sys.path.append(filepath.parent.parent)
# # sys.path.append(filepath.parent.parent.parent)
# sys.path.i
# print(sys.path)


from cascade.xdepth_conversion import XdepthConversion


xh_conv = XdepthConversion()
xh_conv.set_theta(60)
max_xdepth = xh_conv.convert_h2x(0)
# print(xh_conv.convert_h2x(0))



