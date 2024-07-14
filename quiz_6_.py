import pandas as pd
import json

# Function to load JSON data from a file
def load_json_from_file(file_path):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data

# Function to convert the loaded JSON data to a DataFrame
def convert_json_to_df(json_data):
    # Create a list to collect rows
    rows_list = []

    # Process each hit in the results
    for hit in json_data['funfam_scan']['results']:
        # Each hit may contain multiple HSPs (high-scoring pairs), we need to flatten these
        for hsp in hit['hits']:
            # Basic information
            row_data = {
                'query_name': hit['query_name'],
                'match_name': hsp['match_name'],
                'match_description': hsp['match_description'],
                'significance': hsp['significance'],
                'match_type': hsp['match_type'],
                'match_length': hsp['match_length'],
                'match_famnum': hsp['match_funfam_number'],
                'match_taxid': hsp['match_cath_id']['id'],
                'match_model': hsp['hsps'][0]['algorithm'], # Taking algorithm from the first HSP
                'rank': hsp['rank'],
                'data_funam_members': hsp['data']['funfam_members'],
                'term': hsp['data']['rep_id'], # Assuming 'term' correlates to 'rep_id'
                'database_id': hsp['data']['rep_id'],
                'source_id': hsp['data']['rep_source_id']
            }

            # Additional HSP-specific information could be added here

            # Append this row to our list
            rows_list.append(row_data)

    # Create a DataFrame from the list of rows
    df = pd.DataFrame(rows_list)

    return df

# Function to save DataFrame to an Excel file
def save_df_to_excel(df, file_name):
    # Save the DataFrame to an Excel file
    df.to_excel(file_name, index=False)

# Load the JSON data
json_file_path = 'cath_results.JSON' # Adjust the path if the file is in a different directory
json_data = load_json_from_file(json_file_path)

# Convert the JSON to a DataFrame
df = convert_json_to_df(json_data)

# Specify the file name for the Excel file
file_name = 'cath_results.xlsx'

# Save the DataFrame to an Excel file
save_df_to_excel(df, file_name)

# Output the file path
print(file_name)
