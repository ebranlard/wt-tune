import matplotlib as mpl
from matplotlib import pyplot as plt


def dash(p):
    seq = [13, 5]
    p[0].set_dashes(seq)

def dotdash(p):
#     seq = [2, 4, 7, 4]
    seq = [3, 4, 13, 4]
    p[0].set_dashes(seq)

