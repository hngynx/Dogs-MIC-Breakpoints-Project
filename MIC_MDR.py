#!/usr/bin/env python
# coding: utf-8

# # MDR Convert
# 1. If the case is "R" to at least three different classes, then it is marked as "MDR"
# 2. Mark each case into different classes (in bold)
# 3. Add up to see if the sum is up at least 3, mark with "MDR" column
# 
# ## Aminoglycosides
# - Amikacin
# - Gentamicin
# - Neomycin
# 
# ## Beta-Lactams
# - Ampicillin
# - Amoxicillin/Clavulanate (Amoxi/Clav)
# - Piperacillin/Tazobactam (Pipera/Tazob
# - Oxacillin + 2% NaCl (Oxa + 2% NaCl)
# - Cefazolin (1st generation)
# - Cefovecin (3rd generation)
# - Cefpodoxime (3rd generation)
# - Ceftazidime (3rd generation)
# 
# ## Fluoroquinolones
# - Enrofloxacin
# - Marbofloxacin
# - Orbifloxacin
# - Pradofloxacin
# 
# ## Macrolides and Lincosamides
# - Erythromycin (Macrolide)
# - Clindamycin (Lincosamide)
# 
# ## Tetracyclines
# - Doxycycline
# - Minocycline
# - Tetracycline
#  
# ## Sulfonamides
# - Trimethoprim/Sulfamethoxazole (Trim/Sulmeth)
# 
# ## Nitrofurans
# - Nitrofurantoin
# 
# ## Carbapenems
# - Imipenem
# 
# ## Glycopeptides
# - Vancomycin
# 
# ## Chloramphenicol
# - Chloramphenicol
# 

# In[109]:


# Import needed data and packages
import pandas as pd
import numpy as np

# Read data
data=pd.read_excel("/Users/centracy/Desktop/MIC Project/MIC_data.xlsx")
# print(data)
name=data.columns.tolist()[1:77]
data.shape # Checking the rows and columns


# # Check how many resistance in each of the columns

# In[94]:


# Find index of Ampicillin column 
print(data.columns.get_loc("Ampicillin"))

# Look for many "R" it is in from 77th to 94th columns
# data.iloc[:,77:94].dtypes
print((data.iloc[:, 77:94] == "R").sum())


# In[12]:


# Create a column name
n=["Aminoglycosides","Beta-Lactams","Fluoroquinolones",
   "Macro_Linco","Tetracyclines","Sulfonamides","Nitrofurans",
   "Carbapenems","Glycopeptides","Chloramphenicol","MDR"]

# Create new data frame
new=pd.DataFrame(np.zeros((518, 11)), columns=n)
new.head(5)


# In[13]:


# create a loop follows by the citeria
new["Beta-Lactams"] = data[["Ampicillin", 
                            "Amoxi/Clav",
                            "Cefazolin",
                            "Cefovecin",
                            "Cefpodoxime"]].apply(lambda row: (row == "R").any(), axis=1).astype(int)

new["Aminoglycosides"] = data[["Amikacin",
                              "Gentamicin"]].apply(lambda row: (row == "R").any(), axis=1).astype(int)

new["Chloramphenicol"] = (data["Chloramphenicol"] == "R").astype(int)

new["Tetracyclines"] = data[["Doxycycline",
                            "Minocycline",
                            "Tetracycline"]].apply(lambda row: (row == "R").any(), axis=1).astype(int)

new["Fluoroquinolones"] = data[["Orbifloxacin",
                                "Pradofloxacin"]].apply(lambda row: (row == "R").any(), axis=1).astype(int)

new["Sulfonamides"] = (data["Trim/Sulfa"] == "R").astype(int)

new["Nitrofurans"] = (data["Nitrofurantoin"] == "R").astype(int)

new["Carbapenems"] = (data["Imipenem"] == "R").astype(int)

# Mark with MDR if resistant or not
for i in range(len(new)):
    if new.iloc[i, :-1].sum() >= 3:
        new.loc[i, "MDR"] = 1


            
        


# In[95]:


# Delete the R,S,I marked columns before
data=np.delete(data, np.s_[77:93], axis=1)
df=pd.DataFrame(data)

# combined the new data frame with the data now
final=pd.concat([df, new], axis=1)
print(final.head(5))
final.shape


# In[110]:


final.columns.values[1:77] = name
print(len(final.columns[1:77]))  # Columns being replaced
print(len(name))


# In[15]:


print((new.iloc[:,1:11] == 1).sum())


# In[14]:


print(data["Chloramphenicol"].value_counts())
print(new["Chloramphenicol"].value_counts())
# Check the if they matched


# In[111]:


# Export the new data frame out
import openpyxl
final.to_excel("MIC_MDR.xlsx", index=False)


# In[ ]:




