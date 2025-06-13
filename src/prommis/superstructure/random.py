byproduct_opt_conversions = {
    (3, 1): {"Jarosite": 0.75},
    (3, 2): {"Iron oxide": 1},
    (3, 3): {"Residue": 0.25},
    (3, 4): {"Iron hydroxide": 0.5},
    (3, 5): {"Iron oxide": 1},
    (3, 6): {"Iron oxide": 1},
    (5, 4): {"Iron hydroxide": 0.5},
    (5, 5): {"Iron oxide": 1},
    (5, 7): {}
}

for key1, inner_dict in byproduct_opt_conversions.items():
    print(key1)
    for key2, val in inner_dict.items():
        print(key2)
        if val <= 0:
            raise ValueError('error')

        # print(byproduct_opt_conversions[key1][key2])

