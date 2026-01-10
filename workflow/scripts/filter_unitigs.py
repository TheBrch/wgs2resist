import numpy as np
import os, sys
import yaml

with open(os.path.join("config", "config.yaml"), "r") as f:
    config = yaml.load(f, Loader=yaml.FullLoader)


def get_file(i):
    with open(config["metadata"], "r") as f:
        for line in f:
            fields = line.strip().split("\t")
            if i in fields[1]:
                fields[2]
        return i


unitigs = sys.argv[1]

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
