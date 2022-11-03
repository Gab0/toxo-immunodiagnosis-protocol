import numpy as np
import pandas as pd
import sys
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [17.0, 6.4]

fname = sys.argv[1]
output_fname = sys.argv[2]
data = pd.read_csv(fname)
# data = data.sort_values(
#     by="delta_stability_encom",
#     ascending=False,
#     key=abs
# )

DDG = "$\\Delta\\Delta$G"
DDG = "ΔΔG"
rename = {
        "delta_stability_encom": "ΔEstabilidade (kcal/mol)",
        "ddg_prediction": f"Predição {DDG} (kcal/mol)",
        "mutation": "Mutação",
        "delta_vibrational_entropy": "$\\Delta$Vibração Entrópica",
        "dDDG": f"Predição {DDG} (dDDG)"
}


muts_df = pd.read_csv("../muts.csv", keep_default_na=False)
print(muts_df)


def make_graph(mut_names, labels, ys):
    for y in ys:
        plt.plot(mut_names, y)
    plt.legend(labels=labels)
    plt.savefig("compare_ddg.png")


def mut_status(mut):
    possibilities = ["Artificial", "Natural"]
    for row in muts_df.iloc():
        V = row["Variações"]
        if mut in V:
            for pos in possibilities:
                if pos in row["Nome"]:
                    return pos
    return "Unknown"


def parse_dDDG(fname):
    def parse_line(line):
        try:
            return float(line.split(" ")[-1])
        except:
            return None

    with open(fname, encoding="utf8") as F:
        return list(filter(None, map(parse_line, F.readlines())))


def dDDG_to_dyna(mut):
    [c, f, pos, t] = mut.split(" ")
    return "".join([f, pos, t])


def dyna_to_dDDG(mut):
    f = mut[0]
    t = mut[1]
    loc = mut[1:-1]
    return " ".join(["F", f, loc, t])


data["Origem da Mutação"] = data["mutation"].apply(mut_status)

data["dDDG"] = parse_dDDG("dDDG.txt")


prediction_labels = [
    "dDDG",
    "ddg_prediction",
    "delta_stability_encom",
    "delta_vibrational_entropy"
]

prediction_lines = [
    data[label] for label in prediction_labels
]

print(np.cov(prediction_lines))
make_graph(data["mutation"], prediction_labels, prediction_lines)

# Reorder columns.
data = data.sort_values(by="ddg_prediction", ascending=False)
data = data[[
    'mutation',
    'Origem da Mutação',
    'delta_stability_encom',
    'ddg_prediction',
    #'dDDG'
]]

data = data.rename(columns=rename)
print(" & ".join(data.columns))
data.to_csv(output_fname, index=False)
