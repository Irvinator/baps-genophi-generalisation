import pandas as pd
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("--best_hits", required=True)
    parser.add_argument("--out_csv", required=True)

    args = parser.parse_args()

    df = pd.read_csv(args.best_hits, sep="\t", header=None)
    df.columns = ["protein", "cluster"]

    # extract phage id from protein name
    df["phage"] = df["cluster"].str.split("_").str[0]

    # unique proteins per phage
    mapped = df.groupby("phage")["protein"].nunique().reset_index()
    mapped.columns = ["phage", "mapped_proteins"]

    mapped.to_csv(args.out_csv, index=False)

    print("Saved:", args.out_csv)


if __name__ == "__main__":
    main()
