import pandas as pd

# First, we load the data into a pandas DataFrame
df = pd.read_excel('cath_results.xlsx')

# Then we convert 'data_funam_members' to numeric, as it's necessary for the comparison
df['data_funam_members'] = pd.to_numeric(df['data_funam_members'], errors='coerce')

# Now we filter the DataFrame based on the given conditions
df_rna_binding = df[
    (df['match_description'].str.contains('RNA-binding', case=False, na=False)) &
    (df['data_funam_members'] >= 20)
]

# Save the filtered DataFrame to an Excel file called 'rna_binding.xlsx'
df_rna_binding.to_excel('rna_binding.xlsx', index=False)
