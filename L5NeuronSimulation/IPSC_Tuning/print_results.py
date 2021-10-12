import sys
import numpy as np
import pandas as pd

def main():
    ipscs_filename = sys.argv[-1]
    if __file__ != sys.argv[-1]:
        ipscs_filename = sys.argv[-1]
    else:
        raise Exception("Must be run as python print_results <filename>")

    ipsc_csv = pd.read_csv(ipscs_filename)
    ipscs = ipsc_csv["IPSC"]

    print("Mean:", np.mean(ipscs) * 1000, "pA")
    print("Std:", np.std(ipscs) * 1000, "pA")



if __name__ == '__main__':
    main()
