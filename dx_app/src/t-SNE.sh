#!/bin/bash
# t-SNE 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://wiki.dnanexus.com/Developer-Portal for tutorials on how
# to modify this file.

main() {
    #echo "Value of input_sample: '${input_sample[@]}'"
    echo "Value of tissue_type: '${tumor_type}'"

    local_data_dir=$HOME/in
    local_reference_dir=$HOME/reference
    local_output_dir=$HOME/out
    container_data_dir=/data
    container_reference_dir=/reference
    container_output_dir=/results
    metadata_file=/stjude/metadata/combined.csv 
    lookup_file=/stjude/metadata/paper_vs_database_diagnosis_v4_normalized_AlexUpdatesV2_LOOKUP_tSNE.csv
    mkdir $local_reference_dir
    mkdir $local_output_dir

    # Fetch input counts data
    echo ""
    echo "=== Setup ==="
    echo "  [*] Downloading input files ..."
    #dx-download-all-inputs --parallel
    # Get the file ids of the reference count data
    ids=""
    for ((i = 0; i < ${#reference_counts[@]}; i++)) 
    do
        id=${reference_counts[$i]}
        #echo "id: $id"
        clip=$(echo $id | jq '.["$dnanexus_link"]' | sed s'/"//g')
        ids="$ids $clip"
    done
 
    # Setup the download location for reference counts and create the download commands
    # Download with GNU parallel
    mkdir -p $HOME/in/reference_counts/
    echo $ids | xargs -n 100 | sed "s#^#dx download -o $HOME/in/reference_counts/ --no-progress #" > download_all.sh
    parallel --joblog download.log < download_all.sh

    # Look over the inputs, if any, store IDs and download in parallel
    input_ids=""
    if [ ${#input_counts[@]} -gt 0 ]
    then 
       for ((i = 0; i < ${#input_counts[@]}; i++)) 
       do
           id=${input_counts[$i]}
           #echo "id: $id"
           clip=$(echo $id | jq '.["$dnanexus_link"]' | sed s'/"//g')
           input_ids="$input_ids $clip"
       done
    
       mkdir -p $HOME/in/input_counts/
       echo $input_ids | xargs -n 100 | sed "s#^#dx download -o $HOME/in/input_counts/ --no-progress #" > download_inputs.sh
       parallel --joblog download.log < download_inputs.sh
    fi

    echo ""
    echo "  [*] Retrieving covariates for reference data ..."
    covariates_file=$local_data_dir/covariates.txt
    echo -e "Sample\tProtocol\tDiagnosis\tDiagnosisName\tColor" > ${covariates_file}
    
    # Get metadata in parallel for reference data
    echo "Getting metadata for all samples" 
    json=$(echo $ids | xargs python3 /stjude/bin/bulk_describe.py -p $DX_PROJECT_CONTEXT_ID --ids )
    echo $json > metadata.json

    echo "Parsing metadata for each sample"
    echo $json | jq -c '.[] | {name: .name, sample_name: .properties.sample_name, disease: .properties.sj_diseases, type: .properties.sample_type, library: .properties.attr_library_selection_protocol, readlength: .properties.attr_read_length, strandedness: .properties.attr_strandedness, pairing: .properties.attr_read_type, category: .properties.attr_diagnosis_group}' | while read j
    do
      sample_name=$(echo $j | jq -r '.sample_name')

      disease_code=$(echo $j | jq -r '.disease')
      strandedness=$(echo $j | jq -r '.strandedness')
      if [[ "$strandedness" == "Not Available" ]]
      then 
         strandedness="Stranded-Reverse"
      fi 
      if [[ "$strandedness" == "Strandedness unclear" ]]
      then 
         strandedness="Stranded-Reverse"
      fi 
      if [[ "$strandedness" == "undetermined" ]] || [[ "$strandedness" == "Undetermined" ]]
      then 
         echo "Rejecting sample: ${sample_name} [strandedness]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      librarytype=$(echo $j | jq -r '.library')
      readlength=$(echo $j | jq -r '.readlength')
      pairing=$(echo $j | jq -r '.pairing')

      color=$(csvgrep -c 3 -r "^${disease_code}$" $lookup_file |tail -n 1|  csvcut -c 6 -) 
      if [[ "$color" == "Color" ]] || [[ "$color" == "<NA>" ]]
      then
         color='white'
      fi
      disease_name=$(csvgrep -c 3 -r "^${disease_code}$" $lookup_file | tail -n 1 | csvcut -c 2 - | sed 's/"//g') 

      category=$(echo $j | jq -r '.category' )
      if [[ "$category" == "Hematologic Malignancy" ]]
      then
         category="Blood Cancer"
      fi

      protocol="${strandedness}_${librarytype}_${pairing}_${readlength}"

      if [[ "$tumor_type" == "All" ]]  || [[ "$tumor_type" == "$category" ]]
      then
         echo -e "${sample_name}\t${protocol}\t${disease_code}\t${disease_name}\t${color}" | sed 's/"//g' >> ${covariates_file}
      else
         echo "Rejecting sample: ${sample_name} [category]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
      fi
    done

    if [ ${#input_counts[@]} -gt 0 ]
    then 
       echo "Getting metadata for input samples" 
       input_json=$(echo $input_ids | xargs python3 /stjude/bin/bulk_describe.py -p $DX_PROJECT_CONTEXT_ID --ids )
       echo $input_json > input_metadata.json
          
       infile_arg="$container_data_dir/input_counts/*.txt"
       in_arg=""
       echo "Parsing metadata for each input sample"
       set +m
       shopt -s lastpipe
       for ((i = 0; i < ${#input_counts[@]}; i++)) 
       do
       #echo $input_json | jq -c '.[] | {name: .name, sample_name: .properties.sample_name, disease: .properties.sj_diseases, type: .properties.sample_type, library: .properties.library_type, readlength: .properties.read_length, strandedness: .properties.strandedness, pairing: .properties.pairing}' | while read j
       #do
         id=${input_counts[$i]}
         id=$(echo $id | jq -r '.["$dnanexus_link"]')

         echo "id: $id"
         echo "dx describe --json $id"
         j=$(dx describe --json $id)
	 echo "json: $j"
         sample_name=$(echo $j | jq -r '.properties.sample_name')
         disease_code=$(echo $j | jq -r '.properties.disease')
         strandedness=$(echo $j | jq -r '.properties.strandedness')
         if [[ $strandedness == "Not Available" ]]
         then 
            strandedness="Stranded-Reverse"
         fi 
         if [[ $strandedness == "Strandedness unclear" ]]
         then 
            strandedness="Stranded-Reverse"
         fi 
         librarytype=$(echo $j | jq -r '.properties.library')
         readlength=$(echo $j | jq -r '.properties.readlength')
         pairing=$(echo $j | jq -r '.properties.pairing')
         disease_code=$sample_name
         protocol="${strandedness}_${librarytype}_${pairing}_${readlength}"
         echo -e "${sample_name}\t${protocol}\t${disease_code}\t\t" | sed 's/"//g' >> ${covariates_file}
         in_arg="$in_arg --input-sample $sample_name"
         echo "in_arg: $in_arg"
         echo "Adding input sample: $sample_name, $protocol, $disease_code"
       done
       echo ${in_arg}
    fi

    # Fetch gene blacklist
    echo ""
    echo "  [*] Downloading gene blacklist ..."
    dx download -o $local_reference_dir/gene.blacklist.tsv project-F5444K89PZxXjBqVJ3Pp79B4:file-Fk84jFj97xxp1jxP9Zp6JJF4 

    # Fetch Gencode
    echo ""
    echo "  [*] Downloading gencode ..."
    dx download -o $local_reference_dir/gencode.v31.annotation.gtf.gz -r project-F5444K89PZxXjBqVJ3Pp79B4:/pipeline/M2A/gencode.v31.annotation.gtf.gz 
    #dx download -o $local_reference_dir/gencode.v31.annotation.gtf.gz.tbi -r project-F5444K89PZxXjBqVJ3Pp79B4:/pipeline/M2A/gencode.v31.annotation.gtf.gz.tbi 

    # Run interactive t-SNE
    echo ""
    echo "  [*] Running t-SNE ..."
    tissue_arg=
    if [[ "$tumor_type" != "All" ]]
    then
       tissue_arg="--tissue-type \"${tumor_type}\""
    fi
    echo "docker run -v $local_data_dir:$container_data_dir -v $local_reference_dir:$container_reference_dir -v $local_output_dir:$container_output_dir stjudecloud/interactive-tsne:dx_native_app bash -c \"cd $container_output_dir && itsne-main --debug-rscript -b $container_reference_dir/gene.blacklist.tsv -g $container_reference_dir/gencode.v31.annotation.gtf.gz -c $container_data_dir/covariates.txt -o $container_output_dir/${output_name} ${in_arg} ${infile_arg} $container_data_dir/reference_counts/*.txt --save-data ${tissue_arg}\"" 


    docker run -v $local_data_dir:$container_data_dir -v $local_reference_dir:$container_reference_dir -v $local_output_dir:$container_output_dir stjudecloud/interactive-tsne:dx_native_app bash -c "cd $container_output_dir && itsne-main --debug-rscript -b $container_reference_dir/gene.blacklist.tsv -g $container_reference_dir/gencode.v31.annotation.gtf.gz -c $container_data_dir/covariates.txt -o $container_output_dir/${output_name} ${in_arg} ${infile_arg} $container_data_dir/reference_counts/*.txt --save-data ${tissue_arg}" 

    # Upload output  
    tsne_plot=$(dx upload $local_output_dir/${output_name} --brief)
    dx-jobutil-add-output tsne_plot "$tsne_plot" --class=file
    tsne_matrix=$(dx upload $local_output_dir/tsne.txt --brief)
    dx-jobutil-add-output tsne_matrix "$tsne_matrix" --class=file

}
