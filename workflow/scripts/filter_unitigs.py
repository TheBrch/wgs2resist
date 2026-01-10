import numpy as np
import os, sys


unitigs = sys.argv[1]
metadata = sys.argv[2]


def get_file(i):
    i = i.removesuffix("_filtered")
    with open(metadata, "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if i in fields[2]:
                return fields[1]
        return i


with open(unitigs, "r") as f:
    for line in f:
        line = line.rstrip("\n\r")
        parts = line.split("\t")
        if len(parts) == 0:
            continue

        rowname = parts[0]
        numbers = [p for p in parts[1:]]

        try:
            numeric = np.array(numbers, dtype=float)
            n_nonzero = np.count_nonzero(numeric)
            threshold = 0.05 * len(numeric)
            if n_nonzero >= threshold:
                print(line)
        except:
            print("\t".join([get_file(i) for i in parts]))
