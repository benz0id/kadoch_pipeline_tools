"""
BRG1 = '/I'
H3 = '/E'
SESH = "/Users/btudorpr/Desktop/grace struct stuff/7y8r.cxs"
OUT_DIR = "/Users/btudorpr/Desktop/grace struct stuff/7y8r_muts"
d81_pos = '81'
"""

"""
BRG1 = '/A'
H3 = '/K'
SESH = "/Users/btudorpr/Desktop/grace struct stuff/7vdv.cxs"
OUT_DIR = "/Users/btudorpr/Desktop/grace struct stuff/7vdy_muts"
d81_pos = '81'
"""


BRG1 = '/4'
H3 = '/X'
SESH = "/Users/btudorpr/Desktop/grace struct stuff/lab_sesh.cxs"
OUT_DIR = "/Users/btudorpr/Desktop/grace struct stuff/lab_model_muts"
d81_pos = '82'

muts = [
    "K1027E",
    "K1033E",
    "M1036D",
    "M1036K",
    "M1036N",
    "T1034D",
    "T1034V"
]

mut_map = {
"A": "Ala",
"R": "Arg",
"N": "Asn",
"D": "Asp",
"B": "Asx",
"C": "Cys",
"E": "Glu",
"Q": "Gln",
"Z": "Glx",
"G": "Gly",
"H": "His",
"I": "Ile",
"L": "Leu",
"K": "Lys",
"M": "Met",
"F": "Phe",
"P": "Pro",
"S": "Ser",
"T": "Thr",
"W": "Trp",
"Y": "Tyr",
"V": "Val"
}


format_sel = (
f"sel {BRG1}:1027,1033,1034,1035,1036 {H3}:{d81_pos}\n" +
"""
label sel residues text \"{0.label_one_letter_code}{0.number}\" color white bgColor grey height 0.4
addh sel \n""" +
f"sel {BRG1}:1027,1033,1034,1035,1036 {H3}:{d81_pos}\n" +
"""
hide sel atoms
show sel atoms
style sel stick
color sel byhetero
""")

def add_mut(mut: str) -> str:
    i, pos, m = mut[0], mut[1:-1], mut[-1]
    pos = int(pos)
    s = f"""\n
swapaa {BRG1}:{pos} {mut_map[m]} rotLib Dunbrack
color {BRG1}:{pos} byhetero
{format_sel}
sel clear
save "{OUT_DIR}/{mut}.png" supersample 3
open "{SESH}" \n"""
    return s

def add_triple() -> str:
    s = f"""\n
sel
hide sel atoms
label sel delete
sel clear

{format_sel}
sel clear
save "{OUT_DIR}/triple_wt.png" supersample 3

swapaa {BRG1}:1034 ALA rotLib Dunbrack
swapaa {BRG1}:1035 ALA rotLib Dunbrack
swapaa {BRG1}:1036 ALA rotLib Dunbrack
{format_sel}
sel clear
save "{OUT_DIR}/triple.png" supersample 3
open "{SESH}" \n"""
    return s


s = f"""
open "{SESH}"
{format_sel}
sel clear
save "{OUT_DIR}/WT.png" supersample 3
"""

for mut in muts:
    s += add_mut(mut)
s += add_triple()

print(s)