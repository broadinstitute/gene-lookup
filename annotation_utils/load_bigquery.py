"""Load a TSV into BigQuery, overwriting the existing table."""

import argparse
import pandas as pd
from google.cloud import bigquery

PROJECT_ID = "cmg-analysis"
DATASET_ID = "gene_lookup"
TABLE_ID = "combined_gene_disease_association_table"


def main():
    parser = argparse.ArgumentParser(description="Load a TSV into BigQuery.")
    parser.add_argument("tsv_path", help="Path to the TSV file (can be gzipped)")
    args = parser.parse_args()

    print(f"Reading {args.tsv_path}")
    df = pd.read_table(args.tsv_path)
    print(f"Read {len(df)} rows and {len(df.columns)} columns")

    client = bigquery.Client(project=PROJECT_ID)

    # Create dataset if it doesn't exist
    dataset_ref = client.dataset(DATASET_ID)
    try:
        client.get_dataset(dataset_ref)
        print(f"Dataset {DATASET_ID} already exists")
    except Exception:
        dataset = bigquery.Dataset(dataset_ref)
        dataset.location = "US-CENTRAL1"
        client.create_dataset(dataset)
        print(f"Created dataset {DATASET_ID}")

    table_ref = f"{PROJECT_ID}.{DATASET_ID}.{TABLE_ID}"
    print(f"Loading data into {table_ref}")

    job_config = bigquery.LoadJobConfig(write_disposition="WRITE_TRUNCATE")
    job = client.load_table_from_dataframe(df, table_ref, job_config=job_config)
    job.result()

    table = client.get_table(table_ref)
    print(f"Loaded {table.num_rows} rows into {table_ref}")


if __name__ == "__main__":
    main()
