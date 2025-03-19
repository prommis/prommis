### OPEX parameters for iron valorization
# 10 $/kg of jarosite processed
jaro_opex_param = 0.691 * 2.893 * 10

# 5 $/kg of iron hydroxide processed
AFDE_iron_hydrox_opex_param = 0.4740183804 * 1.914 * 5 # AFDE Process
Sel_Leach_iron_hydrox_opex_param = 0.691 * 1.914 * 5 # Selective Leaching Process


print(2.17 + AFDE_iron_hydrox_opex_param)