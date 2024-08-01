# -*- coding: utf-8 -*-
"""
Created on Fri Jun 21 10:26:00 2024

pipeline

@author: Hubert N
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

Hubert Nguyen
"""

import pandas as pd
import os

def excel_formula_pipeline(input_csv, organism_col, start_col, output_csv):
    df = pd.read_csv(input_csv)
    antibiotic_cols = df.columns[df.columns.get_loc(start_col):]

    for col in antibiotic_cols:
        rsi_col_name = col + "_RSI"
        func_name = f"determine_rsi_{col.lower().replace('/', '').replace('-', '')}"

        if func_name in globals():  # Check if function exists
            df[rsi_col_name] = df.apply(
                lambda row: globals()[func_name](row[organism_col], row[col]), axis=1
            )
            insert_position = df.columns.get_loc(col) + 1
            df.insert(insert_position, rsi_col_name, df.pop(rsi_col_name))
        else:
            print(f"No RSI function found for {col}. Skipping.")

    df.to_csv(output_csv, index=False)

# antibiotics for e. coli
# amikacin
def determine_rsi_amikacin(organism, value):
    """
    Determines Amikacin susceptibility (R/I/S) for E. coli based on MIC value.
    """

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
    except ValueError:
        return ""  # Invalid value
    
# ampicillin (UTI)
def determine_rsi_ampicillin(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 8 else "R"
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 8 else "R"  
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff > 8 else "S"
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff > 8 else "S"  
        else:  
            cutoff = float(value)
            return "S" if cutoff <= 8 else "R"

    except ValueError:
        return ""  # Invalid value

# amoxicillin/clavulanate (UTI)
def determine_rsi_amoxiclav(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 8 else "R"
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff < 8 else "R"  # Strictly less than
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff > 8 else "S"
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 8 else "S"
        else:
            cutoff = float(value)
            return "R" if cutoff > 8 else "S"

    except ValueError:
        return ""  # Invalid value

# cefazolin (UTI)
def determine_rsi_cefazolin(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 16 else "R"
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff < 16 else "R"  # Strictly less than
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 32 else "S"
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff > 32 else "S"
        else:
            cutoff = float(value)
            return "R" if cutoff >= 32 else "S"

    except ValueError:
        return ""  # Invalid value

# cefovecin (UTI) 
def determine_rsi_cefovecin(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")

    except ValueError:
        return ""  # Invalid value

# cefpodoxime (UTI)
def determine_rsi_cefpodoxime(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")

    except ValueError:
        return ""  # Invalid value

# ceftazidime (SST)
def determine_rsi_ceftazidime(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")

    except ValueError:
        return ""  # Invalid value

# chloramphenicol
def determine_rsi_chloramphenicol(organism, value):
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""
        
    if value == "<=16":
        return "Invalid"

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 16 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 16 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 32 else ("I" if cutoff == 16 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 32 else ("I" if cutoff == 16 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 32 else ("I" if cutoff == 16 else "S")

    except ValueError:
        return ""
    
# doxycycline
def determine_rsi_doxycycline(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")

    except ValueError:
        return ""  # Invalid value

''' enrofloxacin
def determine_rsi_enrofloxacin(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 0.5 else ("I" if 1 <= cutoff <= 2 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 0.5 else ("I" if 1 <= cutoff <= 2 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 4 else ("I" if 1 <= cutoff <= 2 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 4 else ("I" if 1 <= cutoff <= 2 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 4 else ("I" if 1 <= cutoff <= 2 else "S")

    except ValueError:
        return ""  # Invalid value
    '''

# gentamicin
def determine_rsi_gentamicin(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 2 else ("I" if cutoff == 4 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 8 else ("I" if cutoff == 4 else "S")

    except ValueError:
        return ""  # Invalid value

# imipenem
def determine_rsi_imipenem(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 1 else ("I" if cutoff == 2 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 1 else ("I" if cutoff == 2 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 4 else ("I" if cutoff == 2 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 4 else ("I" if cutoff == 2 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 4 else ("I" if cutoff == 2 else "S")

    except ValueError:
        return ""  # Invalid value

'''# marbofloxacin (SST, UTI)
def determine_rsi_marbofloxacin(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 1 else ("I" if cutoff == 2 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 1 else ("I" if cutoff == 2 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 4 else ("I" if cutoff == 2 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 4 else ("I" if cutoff == 2 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 4 else ("I" if cutoff == 2 else "S")

    except ValueError:
        return ""  # Invalid value
    '''

# minocycline
def determine_rsi_minocycline(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")

    except ValueError:
        return ""  # Invalid value

# nitrofurantoin (UTI)
def determine_rsi_nitrofurantoin(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 32 else ("I" if cutoff == 64 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 32 else ("I" if cutoff == 64 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 128 else ("I" if cutoff == 64 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 128 else ("I" if cutoff == 64 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 128 else ("I" if cutoff == 64 else "S")

    except ValueError:
        return ""  # Invalid value

# orbifloxacin (SST, UTI)
def determine_rsi_orbifloxacin(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 1 else ("I" if 2 <= cutoff <= 4 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 1 else ("I" if 2 <= cutoff <= 4 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 8 else ("I" if 2 <= cutoff <= 4 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 8 else ("I" if 2 <= cutoff <= 4 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 8 else ("I" if 2 <= cutoff <= 4 else "S")

    except ValueError:
        return ""  # Invalid value

# piperacillin/Tazobactam (SST, UTI)
def determine_rsi_pipertazobactam(organism, value):
   
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 8 else ("I" if cutoff == 16 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 8 else ("I" if cutoff == 16 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 32 else ("I" if cutoff == 16 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 32 else ("I" if cutoff == 16 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 32 else ("I" if cutoff == 16 else "S")

    except ValueError:
        return ""  # Invalid value
    
# pradofloxacin (Skin, UTI)
def determine_rsi_pradofloxacin(organism, value):
   
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 0.25 else ("I" if 0.5 <= cutoff <= 1 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 0.25 else ("I" if 0.5 <= cutoff <= 1 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 2 else ("I" if 0.5 <= cutoff <= 1 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 2 else ("I" if 0.5 <= cutoff <= 1 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 2 else ("I" if 0.5 <= cutoff <= 1 else "S")

    except ValueError:
        return ""  # Invalid value
    
#  tetracycline
def determine_rsi_tetracycline(organism, value):
    
    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 4 else ("I" if cutoff == 8 else "R")
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")
        else:
            cutoff = float(value)
            return "R" if cutoff >= 16 else ("I" if cutoff == 8 else "S")

    except ValueError:
        return ""  # Invalid value

# trimethoprim/sulfamethoxazole
def determine_rsi_trimsulfa(organism, value):

    if organism.lower() != "e. coli":
        return ""

    if pd.isna(value) or not isinstance(value, str):
        return ""

    try:
        if "<=" in value:
            cutoff = float(value[2:])
            return "S" if cutoff <= 2 else "R"
        elif "<" in value:
            cutoff = float(value[1:])
            return "S" if cutoff <= 2 else "R"  
        elif ">=" in value:
            cutoff = float(value[2:])
            return "R" if cutoff >= 4 else "S"  
        elif ">" in value:
            cutoff = float(value[1:])
            return "R" if cutoff >= 4 else "S" 
        else:
            cutoff = float(value)
            return "R" if cutoff > 2 else "S" 

    except ValueError:
        return ""  


# input here
input_csv = r"C:\Users\Hubert N\Downloads\test_data.csv" 
organism_col = "Organism Group"
start_col = "Ampicillin" 

# dynamically construct the output file path in the downloads folder
downloads_folder = os.path.expanduser("~\\Downloads")  # downloads folder path
output_filename = "results_rsi.csv"
output_csv = os.path.join(downloads_folder, output_filename)  # full output path

excel_formula_pipeline(input_csv, organism_col, start_col, output_csv)
print(f"Output saved to {output_csv}")

