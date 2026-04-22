"""
Question 2: Data processing
File: simple_somatic_mutation.open.BLCA-CN.tsv.gz
"""
import pandas as pd


def get_allele_mutation_counts(df: pd.DataFrame):
    """
    # Q2.1 – Allele change patterns + unique mutation counts
    Return a table of (mutated_from_allele, mutated_to_allele) patterns
    with the count of unique icgc_mutation_id values for each pattern.
    """
    # Drop duplicate icgc_mutation_id rows first (keep first occurrence)
    unique_mutations = df.drop_duplicates(subset=["icgc_mutation_id"])

    result = (
        unique_mutations
        .groupby(["mutated_from_allele", "mutated_to_allele"])["icgc_mutation_id"]
        .count()
        .reset_index()
        .rename(columns={"icgc_mutation_id": "count_unique_icgc_mutation_id"}) # to match the expected output
        .sort_values(["mutated_from_allele", "mutated_to_allele"])
    )
    return result


def get_sample_mutation_extremes(df: pd.DataFrame):
    """
    # Q2.2 – Sample with highest / lowest unique mutation count
    Find the icgc_sample_id with the highest and lowest count of unique icgc_mutation_id values.
    """
    # Deduplicate: one row per sample, mutation pair
    unique_per_sample = df.drop_duplicates(subset=["icgc_sample_id", "icgc_mutation_id"])

    counts = (
        unique_per_sample
        .groupby("icgc_sample_id")["icgc_mutation_id"]
        .count()
        .reset_index()
        .rename(columns={"icgc_mutation_id": "unique_mutation_count"}) # to match the expected output
    )

    highest = counts.loc[counts["unique_mutation_count"].idxmax()]
    lowest  = counts.loc[counts["unique_mutation_count"].idxmin()]

    return {
        "highest": {"sample_id": highest["icgc_sample_id"],
                    "count":     int(highest["unique_mutation_count"])},
        "lowest":  {"sample_id": lowest["icgc_sample_id"],
                    "count":     int(lowest["unique_mutation_count"])},
    }


if __name__ == "__main__":
    filepath = "COSMIC/simple_somatic_mutation.open.BLCA-CN.tsv.gz"
    # read the file and unzip it
    df = pd.read_csv(filepath, sep="\t", compression="gzip", low_memory=False)

    # embedded print fo a clear answer with the execution 
    print("Q2.1 : Allele change patterns")
    allele_counts = get_allele_mutation_counts(df)
    print(allele_counts.to_string(index=False))

    print("Q2.2 : Sample mutation extremes")
    extremes = get_sample_mutation_extremes(df)
    print(f"Highest unique mutations -> sample: {extremes['highest']['sample_id']}  "
          f"count: {extremes['highest']['count']}")
    print(f"Lowest  unique mutations -> sample: {extremes['lowest']['sample_id']}  "
          f"count: {extremes['lowest']['count']}")