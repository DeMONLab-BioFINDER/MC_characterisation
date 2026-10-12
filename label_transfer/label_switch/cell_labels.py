'''
Label mapping to match SEA-AD PFC dataset to ROSMAP DLPFC that was used before.
'''


celltype_map = {
    # Astrocytes
    "Astro_2": "Ast GRM3",
    "Astro_3": "Ast DPP10",
    "Astro_6-SEAAD": "Ast CHI3L1",
    "Astro_5": "Ast DPP10",
    "Astro_1": "Ast GRM3",
    "Astro_4": "Ast DPP10",

    # OPCs
    "OPC_2": "OPC",
    "OPC_2_2-SEAAD": "OPC",
    "OPC_1": "OPC",

    # Oligodendrocytes
    "Oligo_2": "Oli",
    "Oligo_4": "Oli",
    "Oligo_1": "Oli",
    "Oligo_3": "Oli",

    # Vascular
    "Endo_2": "vasc",
    "VLMC_1": "vasc",
    "Pericyte_1": "vasc",
    "SMC-SEAAD": "vasc",

    # Immune
    "Micro-PVM_2": "Mic P2RY12",
    "Lymphocyte": "T cells",
    "Micro-PVM_3-SEAAD": "Mic TPT1",
    "Micro-PVM_2_1-SEAAD": "Mic MKI67",
    "Micro-PVM_1": "CAMs",
    "Monocyte": "T cells",
    
    #L2/3
    "L2/3 IT_1": "Exc L2-3 CBLN2 LINC02306",
    "L2/3 IT_2": "Exc L2-3 CBLN2 LINC02306",
    "L2/3 IT_3": "Exc L3-4 RORB CUX2",
    "L2/3 IT_5": "Exc L3-4 RORB CUX2",
    "L2/3 IT_6": "Exc L2-3 CBLN2 LINC02306",
    "L2/3 IT_7": "Exc L2-3 CBLN2 LINC02306",
    "L2/3 IT_8": "Exc L2-3 CBLN2 LINC02306",
    "L2/3 IT_10": "Exc L3-4 RORB CUX2", #sus
    "L2/3 IT_12": "Exc L3-4 RORB CUX2",
    "L2/3 IT_13": "Exc L3-4 RORB CUX2",

    #L6
    "L6 IT_1": "Exc L6 THEMIS NFIA",
    "L6 IT_2": "Exc L6 THEMIS NFIA",
    "L6 CT_1": "Exc L6 CT",
    "L6 CT_2": "Exc L6 CT",
    "L6 CT_3": "Exc L6 CT",
    "L6 CT_4": "Exc L6 CT",
    "L6b_1": "Exc L6b",
    "L6b_2": "Exc L6b",
    "L6b_3": "Exc L6b",
    "L6b_4": "Exc L6b",
    "L6b_5": "Exc L6b",
    "L6b_6": "Exc L6b",
    "L6 IT Car3_1": "Exc L5/6 IT Car3'",
    "L6 IT Car3_2": "Exc L5/6 IT Car3'",
    "L6 IT Car3_3": "Exc L5/6 IT Car3'",
    
    #L4
    "L4 IT_1": "Exc L4-5 RORB IL1RAPL2",
    "L4 IT_2": "Exc L4-5 RORB GABRG1",
    "L4 IT_3": "Exc L3-4 RORB CUX2",
    "L4 IT_4": "Exc L3-5 RORB PLCH1",
    
    #L5
    "L5 IT_1": "Exc L5-6 RORB LINC02196",
    "L5 IT_2": "Exc L4-5 RORB IL1RAPL2",
    "L5 IT_3": "Exc L5-6 RORB LINC02196",
    "L5 IT_5": "Exc L4-5 RORB IL1RAPL2",
    "L5 IT_6": "Exc L4-5 RORB IL1RAPL2",
    "L5 IT_7": "Exc L4-5 RORB IL1RAPL2",
    "L5 ET_2": "Exc L5 ET",
    "L5/6 NP_1": "Exc L5/6 NP",
    "L5/6 NP_2": "Exc L5/6 NP",
    "L5/6 NP_3": "Exc L5/6 NP",
    "L5/6 NP_6": "Exc L5/6 NP",
    "L5/6 NP_4": "Exc L5/6 NP",
    
    # Inh LAMP
    "Lamp5_Lhx6_1": "Inh L1-6 LAMP5 CA13",
    "Lamp5_1": "Inh LAMP5 NRG1 (Rosehip)",
    "Lamp5_2": "Inh LAMP5 NRG1 (Rosehip)",
    "Lamp5_3": "Inh LAMP5 NRG1 (Rosehip)",
    "Lamp5_4": "Inh LAMP5 RELN",
    "Lamp5_5": "Inh LAMP5 RELN",
    "Lamp5_6": "Inh LAMP5 RELN",

    # Inh Pax6
    "Pax6_1": "Inh L1 PAX6 CA4",
    "Pax6_2": "Inh L1-2 PAX6 SCGN",
    "Pax6_3": "Inh SGCD PDE3A",
    "Pax6_4": "Inh SGCD PDE3A",

    #Inh Sncg
    "Sncg_1": "Inh PTPRK FAM19A1",
    "Sncg_2": "Inh PTPRK FAM19A1",
    "Sncg_3": "Inh SGCD PDE3A",
    "Sncg_4": "Inh SORCS1 TTN",
    "Sncg_5": "Inh SORCS1 TTN",
    "Sncg_6": "Inh SORCS1 TTN", #???????
    "Sncg_8": "Inh PTPRK FAM19A1",

    #Inh Vip
    "Vip_1": "Inh VIP ABI3BP", # COULD ALSO BE Inh RYR3 TSHZ2
    "Vip_2": "Inh RYR3 TSHZ2", # Could ALSO BE Inh VIP TSHZ2
    "Vip_4": "Inh VIP CLSTN2",
    "Vip_5": "Inh ALCAM TRPM3", # BUT COULD ALSO BE Inh VIP CLSTN2
    "Vip_6": "Inh ALCAM TRPM3",
    "Vip_9": "Inh ALCAM TRPM3",
    "Vip_11": "Inh RYR3 TSHZ2",
    "Vip_12": "Inh VIP CLSTN2", # COULD ALSO BE Inh VIP ABI3BP
    "Vip_13": "Inh RYR3 TSHZ2",
    "Vip_14": "Inh VIP THSD7B",
    "Vip_15": "Inh VIP THSD7B",
    "Vip_16": "Inh VIP THSD7B",
    "Vip_18": "Inh VIP ABI3BP",
    "Vip_19": "Inh VIP ABI3BP",
    "Vip_21": "Inh VIP TSHZ2",
    "Vip_23": "Inh VIP TSHZ2",

    #Inh SST
    "Sst Chodl_1": "Inh L6 SST NPY",
    "Sst Chodl_2": "Inh L6 SST NPY",
    "Sst_1": "Inh L5-6 SST TH",
    "Sst_2": "Inh CUX2 MSR1", # COULD ALSO BE SST MAFB
    "Sst_3": "Inh L3-5 SST MAFB",
    "Sst_4": "Inh FBN2 EPB41L4A",
    "Sst_5": "Inh FBN2 EPB41L4A",
    "Sst_7": "Inh FBN2 EPB41L4A",
    "Sst_9": "Inh ENOX2 SPHKAP",
    "Sst_10": "Inh ENOX2 SPHKAP",
    "Sst_11": "Inh CUX2 MSR1", # COULD ALSO BE SST MAFB
    "Sst_12": "Inh L5-6 SST TH",
    "Sst_13": "Inh ENOX2 SPHKAP", #Could be SST MAFB
    "Sst_19": "Inh L3-5 SST MAFB",
    "Sst_20": "Inh CUX2 MSR1", # COULD ALSO BE SST MAFB
    "Sst_22": "Inh CUX2 MSR1",
    "Sst_23": "Inh CUX2 MSR1",
    "Sst_25": "Inh CUX2 MSR1",

    # Inh PVALB
    "Pvalb_1": "Inh PVALB SULF1",
    "Pvalb_2": "Inh PVALB HTR4",
    "Pvalb_3": "Inh PVALB HTR4",
    "Pvalb_5": "Inh PVALB HTR4",
    "Pvalb_6": "Inh PVALB SULF1",
    "Pvalb_7": "Inh GPC5 RIT2",
    "Pvalb_8": "Inh PVALB SULF1",
    "Pvalb_9": "Inh GPC5 RIT2",
    "Pvalb_10": "Inh L5-6 PVALB STON2",
    "Pvalb_12": "Inh PVALB SULF1",
    "Pvalb_13": "Inh L5-6 PVALB STON2",
    "Pvalb_14": "Inh PVALB SULF1",
    "Pvalb_15": "Inh PVALB HTR4",

    #Inh Chandelier
    "Chandelier_1": "Inh PVALB CA8 (Chandelier)",
    "Chandelier_2": "Inh PVALB CA8 (Chandelier)",
}