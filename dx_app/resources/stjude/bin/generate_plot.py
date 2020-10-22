#!/usr/bin/env python3

import pandas as pd
import json
import argparse
import logzero
import logging

def get_args():
  parser = argparse.ArgumentParser(
      description="Convert t-SNE matrix to HTML plot.")
  parser.add_argument(
      "--tissue-type",
      type=str,
      help="Tumor tissue type: ['Hematologic Malignancy', 'Solid Tumor', 'Brain Tumor']")
  parser.add_argument(
      "--debug",
      default=False,
      action="store_true",
      help="Debug the script.")

  args = parser.parse_args()

  if args.debug:
    logzero.loglevel(logging.DEBUG)
  else:
    logzero.loglevel(logging.INFO)

  return args

def label_samples (row):
  if pd.isna(row['projects']):
      return 'user'
  else:
      return ''

if __name__ == "__main__":
  args = get_args()

  matrix = pd.read_table("tsne.txt")
  matrix = matrix.drop(columns=['classes', 'diagnosisNames', 'color'])
  
  matrix['user'] = matrix.apply(lambda row: label_samples(row), axis=1)
  matrix = matrix.round({'t1': 2, 't2': 2})

  data = matrix.to_csv(sep='\t', index=False, header=False)

  f = open('tsne_pp.html', 'w')

  html = """
<!DOCTYPE html>
<html>
<head>
	 <meta charset="utf-8">
</head>
<body>

<script src="https://proteinpaint.stjude.org/bin/proteinpaint.js" charset="utf-8"></script>

<noscript>This page requires JavaScript. Please enable it in your browser settings.</noscript>

<div id=aaa style="margin:10px"></div>

<script>
runproteinpaint({
	holder:document.getElementById('aaa'),
	parseurl:true,
	nobox:1,
	//block:1,
	genome:'hg38',
	mdssamplescatterplot:{
		dataset:'pediatric-cloud',
		analysisdata:{
            %s
			str:`%s`
        },
		// end of customdata
	}
	/*
	termdb:{
		dev:true,
		state:{
			dslabel:'pediatric-cloud',
			genome:'hg38',
		}
	}
	*/
})
</script>

</body>
</html>
  """
  subset_arg = ''
  if args.tissue_type :
      subset_arg = """			subset:{
				key:'diagnosis_group',
				value:'%s'
			},"""
      subset_arg = subset_arg % (args.tissue_type)
  html = html % (subset_arg, data)
  f.write(html)
  f.close()

