import pandas as pd
from collections import Counter


df = pd.read_csv('CT.matrix.tsv', index_col='Cell.BC')
ct_count = df.sum(axis=1)
ct_count = ct_count.sort_values(axis=0, ascending=False)
print("distribution of celltags per cell:")
print(ct_count.value_counts())
celltags_per_cell = ct_count.value_counts()

mat = pd.read_csv('CT.matrix.tsv', index_col='Cell.BC')
print("nbr of barcodes and celltags:")
print(mat.shape)
print(f'We captured {mat.shape[1]} of 5700 CellTags ' +
      f'in our library ({mat.shape[1]/5700:.2%}).')
mat = mat.loc[(mat != 0).any(axis=1)]
print("distribution of cells per celltag:")
with pd.option_context('display.max_rows', None, 'display.max_columns', None):
    print(mat.sum(axis=0).value_counts())
print(mat.sum(axis=0).value_counts())
cells_per_celltag = mat.sum(axis=0).value_counts()

freqs = {}
for bc in mat.index:
    tmp = mat.loc[[bc]]
    d = tmp.loc[:, (tmp != 0).any()].T.to_dict()
    # d = tmp.loc[:, (tmp != 0).any(axis=0)].T.to_dict()
    freqs[bc] = {k: v for k, v in sorted(d[list(d.keys())[0]].items(),
                                         key=lambda item: item[1], reverse=1)}

# tmp = mat.loc[mat.index[[0]]]
# d = tmp.loc[:, (tmp != 0).any(axis=0)].T.to_dict()
# {k: v for k, v in sorted(d[list(d.keys())[0]].items(), key=lambda item:
    # item[1], reverse=1)}


def id_cases(frequencies):
    cases = {
        'one_CT_freq_1': 0,
        'one_CT_freq_over_1': 0,
        'more_CTs_same_freqs': 0,
        'more_CTs_one_top': 0,
        'more_CTs_approx_same': 0
    }

    for d in list(frequencies.values()):
        freqlist = list(d.values())
        if len(freqlist) == 0:
            print("freqlist shouldn't be empty")
            raise Exception()
        elif len(freqlist) == 1:
            if freqlist[0] == 1:
                cases['one_CT_freq_1'] += 1
            elif freqlist[0] > 1:
                cases['one_CT_freq_over_1'] += 1
        elif len(set(freqlist)) == 1:
            cases['more_CTs_same_freqs'] += 1
        elif freqlist[0] >= 2*freqlist[1] and freqlist[0] - freqlist[1] > 1:
            # maybe some floor stuff or smth here instead ^
            cases['more_CTs_one_top'] += 1
        else:
            cases['more_CTs_approx_same'] += 1
    return cases


cases = id_cases(freqs)
print(cases)

print(f'There are {mat.shape[0]} cells that express at least one CellTag.')


uniques = dict()
mangd = [sorted([y for y in x.keys()]) for x in freqs.values()]
for key in mangd:
    try:
        uniques[tuple(key)] += 1
    except KeyError:
        uniques[tuple(key)] = 1
# uniques
fingerprint_count = dict(Counter(uniques.values()))
fingerprint_count = pd.DataFrame(index=fingerprint_count.keys(),
                                 data={'counts': fingerprint_count.values()})
fingerprint_count.sort_index(inplace=True)
print('Fingerprints based on only CellTags (no frequencies):')
print(fingerprint_count)


uniques = dict()
mangd = [sorted([tuple(y) for y in x.items()]) for x in freqs.values()]
for key in mangd:
    try:
        uniques[tuple(key)] += 1
    except KeyError:
        uniques[tuple(key)] = 1
# uniques
fingerprint_freq_count = dict(Counter(uniques.values()))
fingerprint_freq_count = pd.DataFrame(
    index=fingerprint_freq_count.keys(),
    data={'counts': fingerprint_freq_count.values()})
fingerprint_freq_count.sort_index(inplace=True)
print('Fingerprints based on CellTags with their frequencies:')
print(fingerprint_freq_count)
