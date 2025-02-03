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
  parser.add_argument(
    "--subtype-file",
    type=str,
    help="Subtype definitions file"
  )
  parser.add_argument(
    "--metadata-file",
    type=str,
    help="File containing metadata for the samples"
  )
  parser.add_argument(
    "--tsne-file",
    type=str,
    help="File containing t-SNE matrix for the samples"
  )
  

  args = parser.parse_args()

  if args.debug:
    logzero.loglevel(logging.DEBUG)
  else:
    logzero.loglevel(logging.INFO)

  return args

def label_samples (row):
  if pd.isna(row['projects']) or row['projects'] == 'input':
      return 'user'
  else:
      return None

if __name__ == "__main__":
  args = get_args()

  # Read t-SNE matrix
  matrix = pd.read_table(args.tsne_file)
  # Trim precision on x/y coordinates
  matrix = matrix.round({'t1': 2, 't2': 2})
  # Read the metadata file
  metadata = pd.read_json(args.metadata_file)

  # Convert metadata JSON column 'properties' into a table. 
  # Then combined with the rest of the metadata in a single dataframe.
  prop = pd.json_normalize(metadata['properties'])
  metadata = metadata.drop(columns=['properties'])
  metadata = pd.concat([metadata, prop], axis=1, sort=False)
  combined = matrix.merge(metadata, left_on=['samples','projects'], right_on=['sample_name', 'sj_datasets'], how='left')

  # Get subtype group and label information for each disease code
  subtypes = pd.read_csv(args.subtype_file)
  # Add the subtype information to the dataframe
  combined = combined.merge(subtypes, left_on=['sj_diseases'], right_on=['sj_disease'], how="left")
  
  # Rename t-SNE coordinate columns to x & y. Rename subtype labels to group & label
  combined = combined.rename(columns={"t1":"x", "t2":"y", "t_SNE group": "group","t_SNE_label": "label"})

  combined['highlight'] = combined.apply(lambda row: label_samples(row), axis=1)

  # Drop unneeded columns
  combined = combined.drop(columns=["id", "name", "folder", "sj_ega_accessions", "sj_access_unit", "data_access_level", "sj_dataset_accessions", "vendable", "file_state", "released", "sj_pmid_accessions", "file_type", "sj_embargo_date", "sj_pipeline_name", "sj_pipeline_version", "sj_pub_accessions", "sj_publication_titles", "projects"])
  combined = combined.drop(columns=["sj_disease","sj_long_disease_name_y","diagnosis_group"])

  # Copy user submitted sample names to the standard 'sample_name' column
  combined['sample_name'] = combined['sample_name'].fillna(combined['samples'])

  # Rename duplicated column from merge
  combined = combined.rename(columns={"sj_long_disease_name_x": "sj_long_disease_name"})

  # Drop unneeded columns
  combined = combined.drop(columns=["samples", "classes", "diagnosisNames"])

  # Fill group in a way to avoid duplicates
  combined['group'] = combined['group'].fillna('unknown ' + combined['attr_diagnosis_group'])
  # Fill unknown label values with sj_diseases
  combined['label'] = combined['label'].fillna(combined['sj_diseases'])
  # Fill unknown colors with black
  combined['color'] = combined['color'].replace('Classification_tSNE_color (Manuscript)', 'black')
  # Fill remaining missing metadata values with 'unknown'
  combined = combined.fillna("unknown")

  # Split user samples from reference samples
  user_samples = combined.loc[combined.highlight=='user']
  samples = combined.loc[combined.highlight!='user']

  # Fill any samples that didn't get subtype information with sensible defaults
  samples.fillna({'group':'Other Cancer','label':samples['sj_long_disease_name']}, inplace=True)

  # Create dataframe with just unique diseases for ProteinPaint
  diseases = combined[['sj_diseases','color', 'attr_diagnosis_group', 'group', 'label']]
  diseases = diseases.drop_duplicates()
  diseases['attr_diagnosis_group'] = pd.Categorical(diseases['attr_diagnosis_group'], ["Brain Tumor", "Solid Tumor", "Hematologic Malignancy", "Germ Cell Tumor"])
  diseases = diseases.sort_values(by=['attr_diagnosis_group', 'group', 'label'], key=lambda col: col.str.lower())
  diseases = diseases.set_index('sj_diseases')
  diseases = diseases.dropna()

  # Order group headings [brain, solid, blood]
  groups = combined[['attr_diagnosis_group', 'group']]
  groups = groups.drop_duplicates()
  groups = groups.dropna(how='any', subset={'attr_diagnosis_group', 'group'})
  groups['attr_diagnosis_group'] = pd.Categorical(groups['attr_diagnosis_group'], ["Brain Tumor", "Solid Tumor", "Hematologic Malignancy", "Germ Cell Tumor"])
  groups = groups.sort_values(by=['attr_diagnosis_group', 'group'])
  groups = groups.dropna(how='any', subset={'attr_diagnosis_group', 'group'})
  groups = groups.set_index('group')

  # Create HTML output file and write structure
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
	host:'https://proteinpaint.stjude.org',
	parseurl:true,
	nobox:1,
	//block:1,
	genome:'hg38',
	mdssamplescatterplot:{
		analysisdata:{
            %s, // samplekey
            %s, // attr_levels
            %s, // colorbyattributes
            %s, // sample_attributes
			samples:%s, 
      user_samples:%s
        },
		// end of customdata
	}

})
</script>

</body>
</html>
  """
  
  sample_key = """
  			samplekey:'sample_name'"""
        
  attr_levels = """
			attr_levels:[
				{key:'group'}, // level 1
				{key:'sj_diseases', label: 'label'}, // level 2
			]"""
  
  color = """
			colorbyattributes:[{key:'sj_diseases'}]"""

  sample_attributes = """
			sample_attributes: {
				sj_diseases: {
					label:'Diagnosis',
					values: %s,
        },
        group:{
					label:'Category',
          values: %s
				},
				label:{
					label:'Subtype'
				},
        attr_sex:{
          label: 'Sex'
        },
        attr_age_at_diagnosis:{
          label: 'Age at Diagnosis'
        },
        attr_race:{
          label: 'Race'
        },
        attr_subtype_biomarkers:{
          label: 'Subtype Biomarkers'
        },
        attr_tissue_preservative: {
          label: 'Preservative'
        },
        sj_datasets: {
          label: 'Source'
        },
        sample_type: {
          label: 'Sample Type'
        }
			}"""

  # Fill in disease types in sample attributes
  sample_attributes = sample_attributes % (diseases.to_json(orient="index"), groups.to_json(orient="index"))

  # Fill in additional HTML elements, add reference samples, and user samples
  html = html % (sample_key, attr_levels, color, sample_attributes, samples.to_json(orient="records"), user_samples.to_json(orient="records"))
  f.write(html)
  f.close()

