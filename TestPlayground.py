import pyhrv.time_domain as td
from pyhrv import utils

nn = utils.load_sample_nni()

print(td.sdnn(nn))
print(td.sdsd(nn))