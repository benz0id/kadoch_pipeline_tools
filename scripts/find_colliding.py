seqs = [
    "TATAAT",
    "TCATTC",
    "TCCCGA",
    "TCGAAG",
    "TCGGCA",
    "AGTCAA",
    "CAAAAG",
    "CAACTA",
    "CACCGG",
    "CACTCA",
    "CAGGCG",
    "CATGGC",
    "CATTTT",
    "CACGAT",
    "CCAACA",
    "CGGAAT",
    "CTAGCT",
    "CTATAC",
    "GTGATC",
    "TAATCG",
    "TACAGC",
    "GACGAC"]

def collide(s1, s2, num_diff):
    assert len(s1) == len(s2)
    num_eq = 0
    for i in range(len(s1)):
        if s1[i] == s2[i]:
            num_eq += 1
    return num_eq >= len(s1) - num_diff


def find_collisions(adapters):
    for i, adapter in enumerate(adapters):
        for j, other_adapter in enumerate(adapters[i + 1:]):
            if collide(adapter, other_adapter, 2):
                print(f"Collision! \n\t {i} : {adapter}\n\t {i + j + 1} : {other_adapter}")

find_collisions(seqs)