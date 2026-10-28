#!/usr/bin/env python3
"""Map GmALMT12/QUAC1 (7W6K) functional residues onto the focused alignment and
report conservation in AYB ALMT4_2 and its legume relatives."""
import json
from Bio import SeqIO
import pandas as pd

aln = {r.id: str(r.seq) for r in SeqIO.parse('focused_aln.faa', 'fasta')}
resnums = json.load(open('gmalmt12_resnums.json'))  # seq-pos -> struct resnum
gm = aln['GmALMT12_7W6K_A']
gm_ungapped = gm.replace('-', '')
# struct resnum -> ungapped seq index
rn2idx = {rn: i for i, rn in enumerate(resnums)}

# ungapped index -> alignment column for GmALMT12
idx2col = []
k = 0
for col, c in enumerate(gm):
    if c != '-':
        idx2col.append(col)
        k += 1

def ayb_number(aln_seq, col):
    """1-based residue number in an ungapped sequence at alignment column col."""
    if aln_seq[col] == '-':
        return None
    return sum(1 for c in aln_seq[:col + 1] if c != '-')

# functional residues: (struct_resnum, expected_aa, role)
FUNC = [
    (90, 'W', 'Pore-lining aromatic / gating toggle switch'),
    (100, 'E', 'Pore interhelix network (acidic)'),
    (109, 'K', 'Intramembrane positive ring'),
    (113, 'R', 'Intramembrane positive ring'),
    (158, 'R', 'Intramembrane positive ring'),
    (164, 'K', 'Bottom-of-pore constriction'),
    (165, 'K', 'Bottom-of-pore constriction'),
    (168, 'D', 'Pore interhelix network (acidic)'),
    (169, 'Y', 'Pore interhelix network (aromatic)'),
    (187, 'R', 'Extracellular positive selectivity ring'),
    (198, 'R', 'Extracellular positive selectivity ring'),
    (218, 'P', 'PxWxG layer-interface motif'),
    (220, 'W', 'PxWxG layer-interface motif'),
    (222, 'G', 'PxWxG layer-interface motif'),
    (288, 'W', 'WEP layer-interface motif'),
    (289, 'E', 'WEP layer-interface motif'),
    (290, 'P', 'WEP layer-interface motif'),
]

legumes = ['V7BAP2', 'A0A0R0EAN5', 'I1LEY0', 'A0A4D6N7X3', 'A0A1S2XZY2']
rows = []
for rn, exp, role in FUNC:
    if rn not in rn2idx:
        continue
    col = idx2col[rn2idx[rn]]
    gm_aa = gm[col]
    ayb_aa = aln['AYB_ALMT4_2'][col]
    ayb_num = ayb_number(aln['AYB_ALMT4_2'], col)
    leg_aas = {sp: aln[sp][col] for sp in legumes}
    at9_aa = aln['AtALMT9_Q9LS46'][col]
    # conservation: identical to GmALMT12 residue across AYB+legumes?
    allset = [ayb_aa] + list(leg_aas.values())
    ident = all(a == gm_aa for a in allset)
    # chemically similar classes
    cls = {'K':'+','R':'+','D':'-','E':'-','W':'arom','Y':'arom','F':'arom','P':'P','G':'G'}
    same_class = all(cls.get(a,'o') == cls.get(gm_aa,'x') for a in allset)
    rows.append({
        'GmALMT12_res': f'{gm_aa}{rn}',
        'role': role,
        'AYB_ALMT4_2_res': f'{ayb_aa}{ayb_num}' if ayb_num else '-/gap',
        'AYB_identical': ayb_aa == gm_aa,
        'legumes': ''.join(leg_aas.values()),
        'AtALMT9': at9_aa,
        'identical_all_legume_set': ident,
        'same_physchem_class': same_class,
    })

df = pd.DataFrame(rows)
pd.set_option('display.width', 200); pd.set_option('display.max_columns', 20)
print(df.to_string(index=False))
n_id = df['identical_all_legume_set'].sum()
n_cls = df['same_physchem_class'].sum()
print(f"\nSUMMARY: {n_id}/{len(df)} positions identical across AYB+5 legumes to GmALMT12; "
      f"{n_cls}/{len(df)} conserved by physicochemical class.")
print(f"AYB identical to GmALMT12 at {df['AYB_identical'].sum()}/{len(df)} positions.")
df.to_csv('conservation_table.csv', index=False)
