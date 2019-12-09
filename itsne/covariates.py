import os
import pandas as pd

REQUIRED_HEADERS = ["Sample", "Diagnosis", "Strandedness", "LibraryType", "ReadLength"]

class CovariatesFile:

    def __init__(self, filename):
        self.filename = filename
        self.df = pd.read_csv(self.filename, sep="\t")
        missing_headers = []
        for r in REQUIRED_HEADERS:
            if r not in self.df.columns:
                missing_headers.append(r)

        if len(missing_headers) > 0:
            raise RuntimeError("Covariates file is missing headers: {}. Please see the documentation.".format(", ".join(missing_headers)))
        self.df.set_index('Sample', inplace=True)
        self.df["Covariates"] = self.df["Strandedness"] + "_" + self.df["LibraryType"] + "_" + self.df["ReadLength"]
        del self.df["Strandedness"], self.df["LibraryType"], self.df["ReadLength"]
        self.df = self.df.T