import pandas as pd

reports = [
"SRR11412973_report.txt",
"SRR11412976_report.txt",
"SRR11412979_report.txt",
"SRR11412980_report.txt",
"SRR11412984_report.txt",
"SRR21907296_report.txt",
"SRR21907303_report.txt",
"SRR21907307_report.txt",
"SRR21907330_report.txt",
"SRR21907332_report.txt"
]
column_names = [
    "percentage", 
    "cumulative_count", 
    "direct_count", 
    "tax_rank", 
    "taxid", 
    "taxon"
]

taxonomy = ["G", "F", "O", "C", "P", "D"]
lineage = {rank: None for rank in taxonomy}

for i in range(len(reports)):
    df = pd.read_csv(reports[i], sep="\t", header=None, names=column_names)
    # Ensure the 'percentage' column is numeric
    df["percentage"] = pd.to_numeric(df["percentage"], errors="coerce")

    # Filter for rows where tax_rank is 'S' (species level)
    species_df = df[df["tax_rank"] == "S"]

    # Find the row with the highest percentage
    max_species = species_df.loc[species_df["percentage"].idxmax()]
    species_index = max_species.name

    for idx in range(species_index - 1, -1, -1):
        row_rank = df.loc[idx, "tax_rank"]
        if row_rank in taxonomy and lineage[row_rank] is None:
            lineage[row_rank] = df.loc[idx, "taxon"].strip()
            # If we've found all desired ranks, we can stop.
            if all(lineage[r] is not None for r in taxonomy):
                break

    # Step 3: Create a list of the lineage in the specified order.
    lineage_list = [lineage[rank] for rank in taxonomy]


    # Print the taxon name and its percentage
    #print("Taxon with rank S and highest percentage:")
    print(f"{reports[i].split('_')[0]}: {max_species['taxon'].strip()} ({max_species['percentage']}%)")
    print("Lineage:", lineage_list)
    print()
