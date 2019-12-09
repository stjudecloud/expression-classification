import os
import pandas as pd

NUMBER_OF_COUNT_FILES_READ = 0

class CountFile:

    def __init__(self, filename, verbose=False):
        global NUMBER_OF_COUNT_FILES_READ
        NUMBER_OF_COUNT_FILES_READ += 1
        if verbose:
            print("({}) {}".format(NUMBER_OF_COUNT_FILES_READ, filename), end="\r", flush=True)
        self.filename = filename
        self.df = pd.read_csv(self.filename, sep="\t", header=None)
        name = os.path.basename(self.filename)

        if "." in name:
            name = name.split(".")[0]

        self.df.columns = ['Gene', name]
        self.df.set_index('Gene', inplace=True)
        