import statistics as st

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


def blocking(meta):
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

    patient_list.sort(key=lambda x: (-len(x["samples"]), -x["pos"]))

    positive_patients = sum(x["pos"] > 0 for x in patient_list)
    print(f"Bin count: {positive_patients}")

    bins = [{"size": 0, "pos": 0, "patients": []} for _ in range(positive_patients)]

    for patient in patient_list:
        bins.sort(key=lambda x: (x["pos"], x["size"]))
        bins[0]["size"] += len(patient["samples"])
        bins[0]["pos"] += patient["pos"]
        bins[0]["patients"].append(patient["pid"])

    bin_sizes = [x["size"] for x in bins]
    bin_ratios = [x["pos"] / x["size"] for x in bins]
    bin_stats = {
        "r_mean": st.mean(bin_ratios),
        "r_stdev": st.stdev(bin_ratios),
        "s_mean": st.mean(bin_sizes),
        "s_stdev": st.stdev(bin_sizes),
    }
    bin_stats["s_percent"] = 100 * bin_stats["s_stdev"] / bin_stats["s_mean"]

    return bins, bin_stats


b, b_s = blocking(meta)
print(*b, sep="\n")
print(b_s)
