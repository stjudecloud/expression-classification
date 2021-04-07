#!/usr/bin/env python3

import argparse

import pandas as pd
import numpy as np
from sklearn.neighbors import KDTree

def compute_neighbors(tsne_file):
    # Parse t-SNE coordinates file
    df1 = pd.read_csv(tsne_file, sep='\t')

    # Build a KDTree with the coordinates
    tree = KDTree(df1[['t1', 't2']].values)

    # Pull user input samples to query
    inputs = df1[df1['projects'] == 'input']
    inputs = inputs.reset_index()
    inputs_pos = inputs[['t1', 't2']].values

    summary = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format("sample", "top disease code from n5", "support for top disease code from n5", 'n5 neighbors', 'n5 disease counts', "top disease code from n10", "support for top disease code from n10", 'n10 neighbors', 'n10 disease counts')
    # For each input sample, compute 5 nearest neighbors
    for index, row in inputs.iterrows():
        query = np.reshape(row[['t1','t2']].values, (1,2))
        
        dist, ind = tree.query(query, k=6)
        
        disease = {}
        neighbors = ""
        
        for val in ind[0,:]:
            if not df1.loc[val, 'samples'] == row['samples']: 
                if df1.loc[val, 'classes'] in disease:
                    disease[df1.loc[val, 'classes']] += 1
                else:
                    disease[df1.loc[val, 'classes']] = 1
                if neighbors == "":
                    neighbors += df1.loc[val, 'samples'] + "," + df1.loc[val, 'classes']
                else:
                    neighbors += ";" + df1.loc[val, 'samples'] + "," + df1.loc[val, 'classes']
                    
        counts=""
        top=""
        top_count=""
        for key, value in sorted(disease.items(), key=lambda item: item[1], reverse=True): 
            if counts == "":
                top=key
                top_count=value
                counts += "{0}:{1}".format(key, value)
            else:
                counts += ";{0}:{1}".format(key, value)
            
        dist_n10, ind_n10 = tree.query(query, k=11)
        disease_n10 = {}
        neighbors_n10 = ""

        for val in ind_n10[0,:]:
            if not df1.loc[val, 'samples'] == row['samples']: 
                if df1.loc[val, 'classes'] in disease_n10:
                    disease_n10[df1.loc[val, 'classes']] += 1
                else:
                    disease_n10[df1.loc[val, 'classes']] = 1
                if neighbors_n10 == "":
                    neighbors_n10 += df1.loc[val, 'samples'] + "," + df1.loc[val, 'classes']
                else:
                    neighbors_n10 += ";" + df1.loc[val, 'samples'] + "," + df1.loc[val, 'classes']
                    
        counts_n10=""
        top_n10=""
        top_count_n10=""
        for key, value in sorted(disease_n10.items(), key=lambda item: item[1], reverse=True): 
            if counts_n10 == "":
                top_n10=key
                top_count_n10=value
                counts_n10 += "{0}:{1}".format(key, value)
            else:
                counts_n10 += ";{0}:{1}".format(key, value)
            
        summary += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\n".format(row['samples'], top, top_count, neighbors, counts, top_n10, top_count_n10, neighbors_n10, counts_n10)
        
    return summary

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("tsne", help="t-SNE coordinates file to process")
    parser.add_argument("output_file", default="neighbors.tsv", help="Output file to write neighborhood information")
    args = parser.parse_args()

    neighbors = compute_neighbors(args.tsne)

    with open(args.output_file, "w") as file:
        file.write(neighbors)


if __name__ == "__main__":
    main()