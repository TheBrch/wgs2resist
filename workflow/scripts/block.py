bincount = 5

meta = {
    "s00": {"patient": "P001", "susceptible": 1},
    "s01": {"patient": "P001", "susceptible": 1},
    "s02": {"patient": "P002", "susceptible": 1},
    "s03": {"patient": "P005", "susceptible": 1},
    "s04": {"patient": "P005", "susceptible": 1},
    "s05": {"patient": "P005", "susceptible": 1},
    "s06": {"patient": "P006", "susceptible": 1},
    "s07": {"patient": "P001", "susceptible": 0},
    "s08": {"patient": "P001", "susceptible": 0},
    "s09": {"patient": "P002", "susceptible": 0},
    "s10": {"patient": "P002", "susceptible": 0},
    "s11": {"patient": "P003", "susceptible": 0},
    "s12": {"patient": "P003", "susceptible": 0},
    "s13": {"patient": "P003", "susceptible": 0},
    "s14": {"patient": "P003", "susceptible": 0},
    "s15": {"patient": "P003", "susceptible": 0},
    "s16": {"patient": "P003", "susceptible": 0},
    "s17": {"patient": "P004", "susceptible": 0},
    "s18": {"patient": "P005", "susceptible": 0},
    "s19": {"patient": "P006", "susceptible": 0},
}

# S	p1	p2	p3	p4	p5	p6
# 1	2	1	0	0	3	1
# 0	2	2	6	1	1	1
# 						20

patients = {}
for sample_id, info in meta.items():
    pid = info["patient"]
    if pid not in patients:
        patients[pid] = {"samples": [], "pos": 0, "neg": 0}

    patients[pid]["samples"].append(sample_id)
    if info["susceptible"] == 1:
        patients[pid]["pos"] += 1
    else:
        patients[pid]["neg"] += 1

patient_list = [{**info, "pid": pid} for pid, info in patients.items()]

print(patient_list, "\n")

patient_list.sort(key=lambda x: (-len(x["samples"]), -x["pos"]))

print(patient_list, "\n")

bins = [{"no": i, "size": 0, "pos": 0, "patients": []} for i in range(bincount)]

for patient in patient_list:
    bins.sort(key=lambda x: (x["pos"], x["size"]))
    bins[0]["size"] += len(patient["samples"])
    bins[0]["pos"] += patient["pos"]
    bins[0]["patients"].append(patient["pid"])

print(*bins, sep="\n")
