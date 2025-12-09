#!/usr/bin/env python3


import pandas as pd
import numpy as np
import argparse
import sys

def cluster_by_SNPs(SNPdistance, thresholds=[5, 10, 20, 50, 100],
                    output="clusters", symmetric=True):
    # Load distance matrix (auto-detect delimiter)
    try:
        distMatrix = pd.read_csv(SNPdistance, sep=None, engine='python', index_col=0)
    except Exception as e:
        sys.exit(f"Error reading file {SNPdistance}: {e}")

    # Melt into long format
    mat = distMatrix.copy()
    if symmetric:
        mat.values[np.tril_indices_from(mat, k=0)] = np.nan
        distmat = mat.melt(ignore_index=False, var_name="col", value_name="dist").reset_index().rename(columns={"index": "row"})
        distmat = distmat.dropna(subset=["dist"])
    else:
        distmat = mat.melt(ignore_index=False, var_name="col", value_name="dist").reset_index().rename(columns={"index": "row"})
        distmat = distmat[distmat["row"] != distmat["col"]]
        distmat[["row", "col"]] = np.sort(distmat[["row", "col"]].values, axis=1)
        distmat = distmat.drop_duplicates()

    samples = distMatrix.index.tolist()
    results = pd.DataFrame({"SampleID": samples})

    # Clustering for each threshold
    for thresh in thresholds:
        close_relate = distmat[distmat["dist"] <= thresh]
        names = pd.unique(close_relate[["row", "col"]].values.ravel())

        cluster_results = pd.DataFrame({"SampleID": names,
                                        "Cluster": np.arange(1, len(names)+1)})

        for i in range(len(cluster_results)):
            sample = cluster_results.loc[i, "SampleID"]
            newclose = close_relate[(close_relate["row"] == sample) | (close_relate["col"] == sample)]

            connected = pd.unique(newclose[["row", "col"]].values.ravel())
            clusters = cluster_results.loc[cluster_results["SampleID"].isin(connected), "Cluster"].unique()

            if len(clusters) > 0:
                new_cluster = clusters[0]
                cluster_results.loc[cluster_results["SampleID"].isin(connected), "Cluster"] = new_cluster
                cluster_results.loc[cluster_results["Cluster"].isin(clusters), "Cluster"] = new_cluster

            close_relate = close_relate[~((close_relate["row"] == sample) | (close_relate["col"] == sample))]

        merged = results.merge(cluster_results, on="SampleID", how="left")
        results[f"SNPs{thresh}"] = merged["Cluster"]

    # Rename cluster IDs
    for col in results.columns[1:]:
        unique_clusters = [c for c in results[col].dropna().unique()]
        mapping = {old: f"Cluster_{i+1}" for i, old in enumerate(sorted(unique_clusters))}
        results[col] = results[col].map(mapping)

    # Sort and save
    results = results.sort_values(by=results.columns[1])
    outfile = f"{output}.csv"
    results.to_csv(outfile, index=False)
    print(f"Complete")

    return results


def main():
    parser = argparse.ArgumentParser(description="Cluster samples based on SNP distance matrix.")
    parser.add_argument("--input", required=True, help="Path to SNP distance matrix (CSV or TSV).")
    parser.add_argument("--thresholds", nargs="+", type=float, default=[5,10,20,50,100],
                        help="SNP distance thresholds (space-separated).")
    parser.add_argument("--output", default="clusters", help="Output file prefix.")
    parser.add_argument("--symmetric", action="store_true", help="Indicate if matrix is symmetric.")
    
    args = parser.parse_args()
    cluster_by_SNPs(args.input, thresholds=args.thresholds,
                    output=args.output, symmetric=args.symmetric)


if __name__ == "__main__":
    main()
