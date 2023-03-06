# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
import os
import subprocess
import re
import pandas as pd
import numpy as np
import altair as alt
from altair_saver import save

# %%
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "all"

# %%
project_path = os.getenv("PROJECT_PATH")
assert project_path != None

# %%
cellranger_dir = project_path + "/data/raw/aggregate/2021_083_Qinyu_agg/outs/"
fastq_dir = project_path + "/data/raw/MiSeq_CellTag_Barcode/"
out_dir = project_path + "/data/adhoc/celltag_2023/"
plot_dir = out_dir + "plots/Diverse/"

# %%
# !ls {fastq_dir}

# %%
# !mkdir -p {out_dir}
# !mkdir -p {plot_dir}

# %%
samples = ["Immature", "Diverse"]

# %%
sample = samples[1]
sample

# %%
# !ls {fastq_dir}

# %%
# !zless {fastq_dir}{sample}_S2_L001_R2_001.fastq.gz | head

# %%
# !zless {fastq_dir}{sample}_S2_L001_R2_001.fastq.gz | wc -l

# %%
# !zless {fastq_dir}{sample}_S2_L001_R2_001.fastq.gz | sed "/[^ACTG]/d" | head

# %%
# !zless {fastq_dir}{sample}_S2_L001_R2_001.fastq.gz | sed "/[^ACTG]/d" | wc -l

# %%
# !zless {fastq_dir}{sample}_S2_L001_R2_001.fastq.gz | sed "/[^ACTG]/d" | ag "TGTACG[ACTG]{8}GAATTC" | wc -l

# %%
cmd = f'zless {fastq_dir}{sample}_S2_L001_R2_001.fastq.gz | sed "/[^ACTG]/d" | head'
cmd.split()

# %%
# read_2 = !zless {fastq_dir}{sample}_S2_L001_R{2}_001.fastq.gz | sed "/[^ACTG]/d"

# %%
len(read_2)

# %%
read_2[:10]

# %%
read_2[-10:]

# %%
# read_1 = !zless {fastq_dir}{sample}_S2_L001_R{1}_001.fastq.gz | sed "/[^ACTG]/d"

# %%
len(read_1)

# %%
read_1[:10]

# %%
read_1[-10:]

# %%
# read_1_full = !zless {fastq_dir}{sample}_S2_L001_R{1}_001.fastq.gz
# read_2_full = !zless {fastq_dir}{sample}_S2_L001_R{2}_001.fastq.gz

# %%
print(len(read_1_full))
print(len(read_2_full))

# %%
# read_1 = !zless {fastq_dir}{sample}_S2_L001_R{1}_001.fastq.gz | sed "/[^ACTG]/d"
# read_2 = !zless {fastq_dir}{sample}_S2_L001_R{2}_001.fastq.gz | sed "/[^ACTG]/d"

# %%
print(len(read_1))
print(len(read_2))

# %%
read_1_full[:10]

# %% tags=[]
reads = {}
i = 0
for idx in range(1, len(read_1_full), 4):
    # print(idx)
    reads[i] = {
        "read_1": read_1_full[idx],
        "read_2": read_2_full[idx]
    }
    i += 1
    # print(read_1_full[idx])

# %% tags=[]
reads

# %% tags=[]
for val in reads.values():
    if len(val["read_1"]) < 45 or len(val["read_2"]) < 45:
        print(val)

# %%
len(reads[0]["read_1"])

# %% tags=[]
read_1_for_df = []
read_2_for_df = []
i = 0
for idx in range(1, len(read_1_full), 4):
    read_1_for_df.append(read_1_full[idx])
    read_2_for_df.append(read_2_full[idx])
    i += 1

# %%
df = pd.DataFrame({"read_1": read_1_for_df, "read_2": read_2_for_df})

# %%
df

# %%
df.read_1.apply(len).value_counts()

# %%
df.read_2.apply(len).value_counts()

# %%
df.read_2.apply(
    lambda x:
        len(x)
).value_counts()

# %%
df.read_2[0]

# %%
tmp = re.search("[ACTG]{8}", df.read_2[0])
assert tmp is not None
tmp

# %%
prefix = "TGTACG"
suffix = "GAATTC"

# %%
tmp = re.search(prefix + "([ACTG]{8})" + suffix, df.read_2[0])
assert tmp is not None
tmp

# %%
[x for x in tmp.groups()]

# %%
tmp = re.findall(prefix + "([ACTG]{8})" + suffix, df.read_2[0])
assert tmp is not None
tmp

# %%
df.read_2.apply(
    lambda x:
        re.findall(prefix + "([ACTG]{8})" + suffix, x)
).value_counts()

# %%
df.read_2.apply(
    lambda x:
        len(re.findall(prefix + "([ACTG]{8})" + suffix, x))
).value_counts()

# %%
# df["celltag"] = [x[0] if len(x) > 0 else "" for x in df.read_2.apply(
df["celltag"] = [x[0] if len(x) > 0 else "" for x in df.read_2.apply(
    lambda x:
        re.findall(prefix + "([ACTG]{8})" + suffix, x)
)]

# %%
df

# %%
# !ls {cellranger_dir}

# %%
with open(cellranger_dir + "filtered_feature_bc_matrix/barcodes.tsv", "r") as f:
    barcodes = f.read()[:-1].split("\n")

# %%
barcodes[-5:]

# %%
barcodes = [bc[:-2] for bc in barcodes if bc.endswith("1")]

# %%
barcodes[-5:], len(barcodes)

# %%
# r1_prefix = "TACACGACGCTCTTCCGATCT"
r1_prefix = "GATCT"
r1_suffix = "[T]{2,10}"

# %% tags=[]
sum([len(y) > 0 for y in df.read_1.apply(
    lambda x:
        re.findall(r1_prefix, x)
)])

# %% [markdown]
# How many reads does not have the "Read 1" tag?

# %%
df.shape[0] - sum([len(y) > 0 for y in df.read_1.apply(
    lambda x:
        re.findall(r1_prefix, x)
)])

# %% [markdown]
# ---

# %%
len(barcodes[0])

# %%
df.read_1.apply(
    lambda x:
        len(re.findall(r1_prefix + "([ACTG]{16})", x))
).value_counts()

# %%
df.read_1.apply(
    lambda x:
        len(re.findall(r1_prefix + "([ACTG]{16})" + "([ACTG]{12})" + r1_suffix, x))
).value_counts()

# %%
df.read_1.apply(
    lambda x:
        len(re.findall(r1_prefix + "([ACTG]{16})" + "([ACTG]{12})", x))
).value_counts()

# %%
res_r1 = df.read_1.apply(
    lambda x: 
        re.findall(r1_prefix + "([ACTG]{16})" + "([ACTG]{12})", x)
)
res_r1

# %% tags=[]
res_r1 = [x[0] if len(x) > 0 else "" for x in df.read_1.apply(
    lambda x: 
        re.findall(r1_prefix + "([ACTG]{16})" + "([ACTG]{12})", x)
)]
res_r1

# %%
from collections import Counter

# %%
Counter([len(x) for x in res_r1])

# %%
bcs = []
umis = []
for val in res_r1:
    if len(val) > 0:
        bcs.append(val[0])
        umis.append(val[1])
    else:
        bcs.append("")
        umis.append("")

# %%
df["barcode"] = bcs
df["UMI"] = umis
df

# %%
df.barcode.value_counts()

# %% [markdown]
# df["10X_barcodes"] = df.read_1.apply(
#     lambda x: 
#         re.findall(r1_prefix + "([ACTG]{16})", x)
# )

# %% [markdown]
# df["10X_barcodes"] = [x[0] if len(x) > 0 else "" for x in df.read_1.apply(
#     lambda x: 
#         re.findall(r1_prefix + "([ACTG]{16})", x)
# )]

# %%
# df["10X_barcodes"] = df.read_1.apply(
#     lambda x: (
#         res := re.findall(r1_prefix + "([ACTG]{16})", x),
#         res[0] if len(res) > 0 else ""
#     )
# )

# %% [markdown]
# df["10X_barcodes"]

# %%
df

# %%
#df["10X_barcodes"].value_counts()

# %% [markdown]
# del df["10X_barcodes"]

# %%
df

# %% [markdown]
# ---
#
# # Has a reasonably cleaned dataframe of reads and extracted values

# %% [markdown]
# ## Check 10X barcodes

# %%
Counter([x in barcodes for x in df.barcode])

# %%
Counter([len(x) > 0 for x in df.barcode])

# %%
df["valid_barcode"] = [x in barcodes for x in df.barcode]

# %%
df

# %%
df.valid_barcode.value_counts()

# %% [markdown]
# ## Check celltags

# %%
ref_celltag = project_path + "/data/raw/barcodes_reverse_complementary.csv"
ref_celltag = open(ref_celltag, "r").read().split("\n")[:-1]
ref_celltag = [y[::-1] for y in ref_celltag]
ref_celltag[:8], len(ref_celltag)


# %%
def nucleobase_complement(base):
    if base == "T":
        return "A"
    elif base == "A":
        return "T"
    elif base == "C":
        return "G"
    elif base == "G":
        return "C"

def compliment(sequence):
    return "".join([nucleobase_complement(base) for base in sequence])

assert compliment("ACTG") == "TGAC"
    
ref_celltag[0]
compliment(ref_celltag[0])

# %%
ref_celltag[:3]
ref_celltag = [compliment(ct) for ct in ref_celltag]
ref_celltag[:3]

# %%
Counter([x in ref_celltag for x in df.celltag])

# %%
df["valid_celltag"] = [x in ref_celltag for x in df.celltag]

# %%
df

# %%
df[["valid_barcode", "valid_celltag"]]

# %%
(df["valid_barcode"] & df["valid_celltag"]).value_counts()

# %% [markdown]
# # Plots

# %% [markdown]
# ### Valid barcodes _and_ celltags

# %%
source = (df["valid_barcode"] & df["valid_celltag"]).value_counts()
source = list(source)
source = pd.DataFrame(data = {"values": source, "Classification": ["Invalid barcode + celltag", "Valid barcode + celltag"]})
source

# %%
base = alt.Chart(source).encode(
    theta=alt.Theta(field="values", type="quantitative", stack=True),
    color=alt.Color(field="Classification", type="nominal"),
)
pie = base.mark_arc(outerRadius=120)
text = base.mark_text(radius=140, size = 15).encode(text="values:Q")

plot_bcct = pie + text
plot_bcct = plot_bcct.properties(
    title=[
        "Number of reads with CellTags+10X barcodes",
        "that exist in provided CellTag whitelist",
        "and CellRanger filtered barcodes.tsv"]
)

# %% [markdown]
# ### Valid barcodes

# %%
source = df["valid_barcode"].value_counts()
source = list(source)
source = pd.DataFrame(data = {"values": source, "Classification": ["Invalid barcode", "Valid barcode"]})
source

# %%
base = alt.Chart(source).encode(
    theta=alt.Theta(field="values", type="quantitative", stack=True),
    color=alt.Color(field="Classification", type="nominal"),
)
pie = base.mark_arc(outerRadius=120)
text = base.mark_text(radius=150, size = 15).encode(text="values:Q")

plot_bc = pie + text
plot_bc = plot_bc.properties(
    title=["Number of reads with 10X barcodes that exist",  "in CellRanger filtered barcodes.tsv"]
)

# %% [markdown]
# ### Valid CellTags

# %%
source = df["valid_celltag"].value_counts()
source = list(source)
source = pd.DataFrame(data = {"values": source, "Classification": ["Invalid celltag", "Valid celltag"]})
source

# %%
base = alt.Chart(source).encode(
    theta=alt.Theta(field="values", type="quantitative", stack=True),
    color=alt.Color(field="Classification", type="nominal"),
)
pie = base.mark_arc(outerRadius=120)
text = base.mark_text(radius=145, size = 15).encode(text="values:Q")

plot_ct = pie + text
plot_ct = plot_ct.properties(
    title=["Number of reads with CellTags", "that exist in provided CellTag whitelist"]
)

# %% [markdown]
# ### Valid barcodes vs. 10X barcodes

# %%
source = (df["valid_barcode"] & df["valid_celltag"]).value_counts()
n_valid = int(list(source)[1])
n_valid

# %%
len(barcodes)

# %%
source = pd.DataFrame(
    data = {"values": [n_valid, len(barcodes)-n_valid],
            "Classification": ["Valid bc+ct", "Invalid bc+ct"]})
source

# %%
base = alt.Chart(source).encode(
    theta=alt.Theta(field="values", type="quantitative", stack=True),
    color=alt.Color(field="Classification", type="nominal"),
)
pie = base.mark_arc(outerRadius=120)
text = base.mark_text(radius=150, size = 15).encode(text="values:Q")

plot_yield = pie + text
plot_yield = plot_yield.properties(
    title="Yield of CellTags, compared to valid barcodes from CellRanger (!contains duplicates atm!)"
)

# %% [markdown]
# ### Concatenate plots

# %%
plot_bc
plot_ct
plot_bcct
# plot_yield

# %%
save(plot_bc, plot_dir + "pie_valid_barcodes.svg")
save(plot_ct, plot_dir + "pie_valid_celltags.svg")
save(plot_bcct, plot_dir + "pie_valid_celltags_and_barcodes.svg")

# %% [markdown]
# ## Check duplicates

# %%
len(barcodes)

# %%
df[df["valid_barcode"] & df["valid_celltag"]]["barcode"].value_counts()

# %%
len(df[df["valid_barcode"] & df["valid_celltag"]]["barcode"].value_counts()) / len(barcodes)

# %%
valids = df[df["valid_barcode"] & df["valid_celltag"]]
valids = valids.reset_index(names="old_index")
valids = valids[valids.columns.tolist()[3:] + valids.columns.tolist()[:3]]
valids

# %%
source = dict(valids.barcode.value_counts())
source = pd.DataFrame(data={"barcode": source.keys(), "values": source.values()})
source = source.sort_values(by="values", ascending=False)
source

# %%
# alt.Chart(source).mark_bar().encode(x="barcode:N", y="values:Q")
with alt.data_transformers.disable_max_rows():
    alt.Chart(source).mark_bar().encode(
        y="count():Q", x=alt.X("values:Q", bin=alt.Bin(maxbins=60))
    )

# %%
">=10:", source[source["values"] >= 10].shape[0], "<10:", source[source["values"] < 10].shape[0]

# %%
source[source["values"] >= 10].shape[0] / len(barcodes)

# %%
f"Number of 10X barcodes in the count matrix: {len(barcodes)}"

# %%
old_source = source

# %%
source = pd.DataFrame(
    data={
        "values": [
            old_source[old_source["values"] >= 10].shape[0],
            old_source[old_source["values"] < 10].shape[0],
        ],
        "Classification": [">= 10", "< 10"],
    }
)
source

# %%
base = alt.Chart(source).encode(
    theta=alt.Theta(field="values", type="quantitative", stack=True),
    color=alt.Color(field="Classification", type="nominal"),
)
pie = base.mark_arc(outerRadius=120)
text = base.mark_text(radius=150, size = 15).encode(text="values:Q")

plot = pie + text
plot = plot.properties(
    title=["Number of read 10X barcodes with a number", "of associated celltags (thresholded)"]
)
plot
save(plot, plot_dir + "nbr_highqual_10Xbarcodes.svg")

# %%
source = pd.DataFrame(
    data={
        "values": [
            old_source[old_source["values"] >= 10].shape[0] / len(barcodes),
            (1 - (old_source[old_source["values"] >= 10].shape[0] / len(barcodes))),
        ],
        "Classification": [["\"High-quality\" read", "10X bc's (>= 10)"], "Total # 10X bc's"],
    }
)
source

# %%
base = alt.Chart(source).encode(
    theta=alt.Theta(field="values", type="quantitative", stack=True),
    color=alt.Color(field="Classification", type="nominal"),
)
pie = base.mark_arc(outerRadius=120)
text = base.mark_text(radius=155, size=15).encode(
    text=alt.Text("values:Q", format=".2%")
)

plot = pie + text
plot = plot.properties(
    title=[
        f"Fraction of read \"high-quality\" 10X barcodes ({old_source[old_source['values'] >= 10].shape[0]})",
        f"of total number of 10X barcodes ({len(barcodes)})",
    ]
)
plot
save(plot, plot_dir + "fraction_of_all_highqual_10Xbarcodes.svg")

# %%
alt.Chart.configure_mark

# %%
counts = dict(valids.barcode.value_counts())

# %%
ct_counts = {}
for x in [x for x in list(dict(valids.barcode.value_counts()).keys())]:
    ct_counts[x] = valids[valids.barcode == x].celltag.value_counts()

# %%
len(ct_counts)

# %% [markdown]
# How many unique celltag reads for each valid barcode?

# %%
Counter([len(x) for x in ct_counts.values()])

# %% tags=[]
ratios = {}
for k, v in {k: v for k, v in ct_counts.items() if len(v) > 1}.items():
    ratios[k] = v[0] / sum(v[1:])
ratios

# %%
source = pd.DataFrame({"barcodes": ratios.keys(), "values": ratios.values()})
source

# %%
from altair import expr, datum
# alt.Chart(source).mark_bar().encode(x="barcodes", y="values")
with alt.data_transformers.disable_max_rows():
    plot = alt.vconcat(
        alt.Chart(source[source["values"] < 6])
        .mark_bar()
        .encode(y="count():Q", x=alt.X("values:Q", bin=alt.Bin(maxbins=60))),
        alt.Chart(source)
        .mark_bar()
        .encode(y="count():Q", x=alt.X("values:Q", bin=alt.Bin(maxbins=90)))
    ).properties(
        title=[
            "Ratio of first celltag count vs the rest,",
            "for barcodes with more than one associated celltag",
        ]
    )
    plot
    save(plot, plot_dir + "ratio_one_first.svg")

# %% [markdown]
# How big a part of all barcodes have more than one celltag?

# %%
source.shape[0] / len(ct_counts.keys())

# %% [markdown]
# How big a part of the barcodes with more than one celltag are hard to define? I.e. ratio 1 between max celltag and all others associated with the barcode.

# %%
source.values[source["values"] <= 1].shape, source.shape, source.values[source["values"] <= 1].shape[0] / source.shape[0]

# %% [markdown]
# How big a part of all barcodes are "hard to define"?

# %%
source.values[source["values"] <= 1].shape[0] / len(ct_counts.keys())

# %%

# %% tags=[]
ratios_two_top = {}
for k, v in {k: v for k, v in ct_counts.items() if len(v) > 1}.items():
    if len(v) > 2:
        ratios_two_top[k] = sum(v[0:1]) / sum(v[2:])
    else:
        ratios_two_top[k] = v[0] / sum(v[1:])
ratios_two_top

# %%
source = pd.DataFrame({"barcodes": ratios_two_top.keys(), "values": ratios_two_top.values()})
source

# %%
from altair import expr, datum
# alt.Chart(source).mark_bar().encode(x="barcodes", y="values")
with alt.data_transformers.disable_max_rows():
    plot = alt.vconcat(
        alt.Chart(source[source["values"] < 6])
        .mark_bar()
        .encode(y="count():Q", x=alt.X("values:Q", bin=alt.Bin(maxbins=60))),
        alt.Chart(source)
        .mark_bar()
        .encode(y="count():Q", x=alt.X("values:Q", bin=alt.Bin(maxbins=90)))
    ).properties(
        title=[
            "Ratio of first two celltags count vs the rest,",
            "for barcodes with more than one associated celltag",
        ]
    )
    plot
    save(plot, plot_dir + "ratio_two_first.svg")

# %%

# %%
for x in [x for x in list(dict(valids.barcode.value_counts()).keys())[:30]]:
    print(">>", x, "<<")
    print(valids[valids.barcode == x].celltag.value_counts())

# %%
