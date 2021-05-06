import random as rd
import pandas as pd
import numpy as np
import sys

# def main2(iseed):
#     np.random.seed(iseed.astype(int))
#     return np.random.rand()


def main(iseed, u):
    np.random.seed(iseed)
    return np.random.rand() + u


if __name__ == '__main__':
    iseed = int(sys.argv[1])
    u = int(sys.argv[2])
    sys.stdout.write(str(main(iseed, u)))
