'''
Cell subtype mapping to be applied.
'''
# Combined mappings for cell subtypes
cell_subtype_mappings = {
    # Vascular cells
    'Fib_SLC4A4': 'Vasc',
    'Fib_FLRT2': 'Vasc',
    'Per': 'Vasc',
    'End': 'Vasc',
    'SMC': 'Vasc',

    # Inhibitory neurons Group 1
    'Inh_PTPRK_FAM19A1': 'Inh_G1', #closest
    'Inh_SORCS1_TTN': 'Inh_G1', #closes
    'Inh_LAMP5_RELN': 'Inh_G1', #less close, still close tho
    
    # Inhibitory neurons Group 2 (actually group 1 now)
    'Inh_L12_PAX6_SCGN': 'Inh_G1',
    'Inh_L1_PAX6_CA4': 'Inh_G1',
    'Inh_SGCD_PDE3A': 'Inh_G1',


    # Inhibitory neurons Group 3
    'Inh_LAMP5_NRG1': 'Inh_G3',  #has 26k cells, one distinct group with this other one that:
    'Inh_L16_LAMP5_CA13': 'Inh_G3', #has 15k cells only


    # Inhibitory neurons Group 4
    'Inh_RYR3_TSHZ2': 'Inh_G4', #24k cells, close
    'Inh_ALCAM_TRPM3': 'Inh_G4', #10k cells, close

    'Inh_VIP_THSD7B': 'Inh_G4', #10k cells, less close
    
    
    # Inhibitory neurons Group 5
    'Inh_VIP_CLSTN2': 'Inh_G5', #17k cells, close
    'Inh_VIP_ABI3BP': 'Inh_G5', #14k cells, close 

    'Inh_VIP_TSHZ2': 'Inh_G5', #12k cells, less close    
    
    
    # Inhibitory neurons Group 6
    'Inh_L56_SST_TH': 'Inh_G6', #4k cells, close
    'Inh_FBN2_EPB41L4A': 'Inh_G6', #6k cells, close 

    'Inh_L35_SST_MAFB': 'Inh_G6', #23k cells, less close    
    'Inh_ENOX2_SPHKAP': 'Inh_G6', #12k cells, less close
    'Inh_L6_SST_NPY' : 'Inh_G6',    
    
    
    # Inhibitory neurons Group 7
    'Inh_L56_PVALB_STON2': 'Inh_G7', #5k cells, close
    'Inh_GPC5_RIT2': 'Inh_G7', #5k cells, close 
    
    'Inh_PVALB_HTR4': 'Inh_G7', #42k cells, less close    
    'Inh_PVALB_SULF1': 'Inh_G7', #21k cells, less close
    
    
    # Inhibitory neurons Group 8
    'Inh_PVALB_CA8': 'Inh_G8', #15k cells, less close    
    'Inh_CUX2_MSR1': 'Inh_G8', #21k cells, less close
    
    
    # Immune cells
    'Mic_TPT1': 'Mic_all',
    'Mic_P2RY12': 'Mic_all',
    'Mic_MKI67': 'Mic_all',
    'T_cells': 'immun',
    'CAMs': 'immun',
    
    # Excitatory neurons Group1
    'Exc_L56_IT_Car3' : 'Exc_G1',
    'Exc_L5_ET' : 'Exc_G1'
}
