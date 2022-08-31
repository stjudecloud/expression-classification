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

set -e -o pipefail

main() {
   echo "Value of tissue_type: '${tumor_type}'"
   echo "Value of all_strandedness: '${all_strandedness}'"
   echo "Value of all_library_type: '${all_library_type}'"
   echo "Value of all_read_length: '${all_read_length}'"
   echo "Value of all_pairing: '${all_pairing}'"
   echo "Value of include_pdx: '${include_pdx}'"
   echo "Value of output_name: '${output_name}'"
   echo "Value of intermediate file: '${intermediate}'"
   echo "Value of gene list: '${gene_list}'"
   echo "Value of preservatives: '${preservatives}'"

   output_prefix=$(basename ${output_name} ".html")

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
   echo "  [*] Check for duplicate files ..."

   # Get file names for reference count files
   file_names=()
   for ((i = 0; i < ${#reference_counts_name[@]}; i++))
   do
      name=${reference_counts_name[$i]}
      file_names+=( $name )
   done
   echo "   Getting unique file count"
   uniqueNum=$(printf '%s\n' "${file_names[@]}"|awk '!($0 in seen){seen[$0];c++} END {print c}')
   duplicateFiles=$(printf '%s\n' "${file_names[@]}"|awk '!($0 in seen){seen[$0];next} 1' | xargs echo)

   (( 0 == ${#file_names[@]} )) && echo "Found no reference files" && \
   echo "{\"error\": {\"type\": \"AppError\", \"message\": \"No reference files found.\"}}" > job_error.json && \
   exit 1


   echo "   Checking unique file count"
   (( $uniqueNum != ${#file_names[@]} )) && echo "Found duplicates" && \
   echo $duplicateFiles && \
   echo "Consider deduplicating with https://github.com/stjudecloud/utilities/blob/master/stjudecloud_utilities/deduplicate_feature_counts.py" && \
   echo "{\"error\": {\"type\": \"AppError\", \"message\": \"Duplicate file names found. $duplicateFiles. Consider deduplicating with https://github.com/stjudecloud/utilities/blob/master/stjudecloud_utilities/deduplicate_feature_counts.py.\"}}" > job_error.json && \
   exit 1

   echo "  [*] Downloading input files ..."
   # Get the file ids of the reference count data
   ids=""
   for ((i = 0; i < ${#reference_counts[@]}; i++))
   do
      id=${reference_counts[$i]}
      clip=$(echo $id | jq '.["$dnanexus_link"]' | sed s'/"//g')
      ids="$ids $clip"
   done

   # Setup the download location for reference counts and create the download commands
   # Download with GNU parallel
   mkdir -p $HOME/in/reference_counts/
   echo $ids | xargs -n 1 | sed "s#^#dx download -f -o $HOME/in/reference_counts/ --no-progress ${DX_PROJECT_CONTEXT_ID}:#" > download_all.sh
   parallel --retries 20 --results download_outputs --joblog download.log < download_all.sh > download.stdout

   # Loop over the inputs, if any, store IDs and download in parallel
   input_ids=""
   if [ ${#input_counts[@]} -gt 0 ]
   then
      for ((i = 0; i < ${#input_counts[@]}; i++))
      do
         id=${input_counts[$i]}
         clip=$(echo $id | jq '.["$dnanexus_link"]' | sed s'/"//g')
         input_ids="$input_ids $clip"
      done

      mkdir -p $HOME/in/input_counts/
      echo $input_ids | xargs -n 100 | sed "s#^#dx download -o $HOME/in/input_counts/ --no-progress #" > download_inputs.sh
      parallel --joblog download.log < download_inputs.sh
   fi

   # Metadata functions
   function get_property(){
      local prop=$(echo $1 | jq -r ".properties.$2")
      echo ${prop}
   }
   function get_sample_name() {
      echo $(get_property "$1" 'sample_name')
   }
   function get_disease_code() {
      echo $(get_property "$1" 'sj_diseases')
   }
   function get_project() {
      echo $(get_property "$1" 'sj_datasets')
   }
   function get_strandedness(){
      echo $(get_property "$1" 'attr_lab_strandedness')
   }
   function get_library(){
      echo $(get_property "$1" 'attr_library_selection_protocol')
   }
   function get_readlen(){
      echo $(get_property "$1" 'attr_read_length')
   }
   function get_pairing(){
      echo $(get_property "$1" 'attr_read_type')
   }
   function get_platform(){
      echo $(get_property "$1" 'attr_sequencing_platform')
   }
   function get_type(){
      echo $(get_property "$1" 'sample_type')
   }
   function get_preservative(){
      echo $(get_property "$1" 'attr_tissue_preservative')
   }
   function get_category(){
      echo $(get_property "$1" 'attr_diagnosis_group')
   }

   echo ""
   echo "  [*] Retrieving covariates for reference data ..."
   covariates_file=$local_data_dir/covariates.txt
   echo -e "Sample\tProtocol\tDiagnosis\tDiagnosisName\tColor\tProjects\tPreservative" > ${covariates_file}

   # Get metadata in parallel for reference data
   echo "Getting metadata for all samples"
   json=$(echo $ids | xargs python3 /stjude/bin/bulk_describe.py -p $DX_PROJECT_CONTEXT_ID --ids )
   echo $json > metadata.json

   # Prepare filters
   excluded_types='germline|cell line'
   if [ "${include_pdx}" != "true" ]
   then
      excluded_types="${excluded_types}|xenograft"
   fi

   # retrieve and set excluded_preservatives varaible (from preservative selected in dxapp.json file
   excluded_preservatives="Not Available"
   if [ "$preservatives" == "Fresh/Frozen" ]
   then
      excluded_preservatives="${excluded_preservatives}|FFPE"
   elif [ "$preservatives" == "FFPE" ]
   then
      excluded_preservatives="${excluded_preservatives}|Fresh/Frozen"
   fi 

   # Parse metadata for reference files
   echo "Parsing metadata for each sample"
   echo $json | jq -c '.[]' | while read j
   do
      sample_name=$(get_sample_name "$j")

      disease_code=$(get_disease_code "$j")

      projects=$(get_project "$j")

      # Strandedness: [Unstranded, Stranded-Forward, Stranded-Reverse]
      strandedness=$(get_strandedness "$j")

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

      # Library type: [total, polyA]
      librarytype=$(get_library "$j")
      if [[ "$librarytype" == "Not Available" ]]
      then
         echo "Rejecting sample: ${sample_name} [library type]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      # Read length in base pairs
      readlength=$(get_readlen "$j")
      if [[ "$readlength" == "Not Available" ]] || [ $(echo $readlength | grep -c "and") -gt 0 ]
      then
         echo "Rejecting sample: ${sample_name} [read length]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      # Read pairing: [Paired-end, Single-end]
      pairing=$(get_pairing "$j")
      if [[ "$pairing" == "Not Available" ]]
      then
         echo "Rejecting sample: ${sample_name} [read pairing]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      # Remove machines: MiSeq, SOLEXA, or unknown
      # 2020-11-09: AG - keep mixed machines
      platform=$(get_platform "$j")
      if [ $(echo $platform | grep -cE 'MiSeq|SOLEXA|Not Available') -gt 0 ]
      then
         echo "Rejecting sample: ${sample_name} [platform]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      # Check for sample type. Remove germline, cell line, xenograft
      type=$(get_type "$j")
      if [ $(echo $type | grep -cE "${excluded_types}") -gt 0 ]
      then
         echo "Rejecting sample: ${sample_name} [sample type]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      # Remove user defined samples with selected preservatives
      preservative=$(get_preservative "$j")
      if [ $(echo $preservative | grep -cE "${excluded_preservatives}") -gt 0 ]
      then
         echo "Rejecting sample: ${sample_name} [preservative]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
         continue
      fi

      # Lookup color code by disease code
      color=$(csvgrep -c 3 -r "^${disease_code}$" $lookup_file |tail -n 1|  csvcut -c 6 -) 
      if [[ "$color" == "Color" ]] || [[ "$color" == "<NA>" ]]
      then
         color='#aba9a9'
      fi
      # Lookup normalized long disease name by disease code
      disease_name=$(csvgrep -c 3 -r "^${disease_code}$" $lookup_file | tail -n 1 | csvcut -c 2 - | sed 's/"//g') 

      category=$(get_category "$j")
      if [[ "$category" == "Hematologic Malignancy" ]]
      then
         category="Blood Cancer"
      fi

      protocol="${strandedness}_${librarytype}_${pairing}_${readlength}"

      if [[ "$tumor_type" == "All" ]]  || [[ "$tumor_type" == "$category" ]]
      then
         echo -e "${sample_name}\t${protocol}\t${disease_code}\t${disease_name}\t${color}\t${projects}\t${preservative}" | sed 's/"//g' >> ${covariates_file}
      else
         echo "Rejecting sample: ${sample_name} [category]"
         file_name=$(echo $j | jq '.name' | sed 's/\"//g')
         rm $HOME/in/reference_counts/${file_name}
      fi
   done

   # Handle Dana Farber PDX samples
   if [ "${include_pdx}" == "true" ]
   then
      # Get the file ids of the reference count data
      dfci_ids=""
      for file in $(dx ls "project-F5444K89PZxXjBqVJ3Pp79B4:/pipeline/RNA-Seq Expression Classification/dana_farber_pdx/")
      do
         json=$(dx describe --json "project-F5444K89PZxXjBqVJ3Pp79B4:/pipeline/RNA-Seq Expression Classification/dana_farber_pdx/${file}")
         id=$(echo $json | jq -r '.id')
         # Download HTSeq count file
         dx download -f -o $HOME/in/reference_counts/ --no-progress project-F5444K89PZxXjBqVJ3Pp79B4:${id}
         dfci_ids="${dfci_ids} ${id}"

         # Get metadata entries
         sample_name=$(get_sample_name "$json")
         disease_code=$(get_disease_code "$json")
         projects=$(get_project "$json")
         strandedness=$(get_strandedness "$json")
         librarytype=$(get_library "$json")
         readlength=$(get_readlen "$json")
         pairing=$(get_pairing "$json")

         # look up color and disease name values
         color=$(csvgrep -c 3 -r "^${disease_code}$" $lookup_file |tail -n 1|  csvcut -c 6 -)
         if [[ "$color" == "Color" ]] || [[ "$color" == "<NA>" ]]
         then
            color='#aba9a9'
         fi
         disease_name=$(csvgrep -c 3 -r "^${disease_code}$" $lookup_file | tail -n 1 | csvcut -c 2 - | sed 's/"//g')

         # build protocol string and write covariates to file
         protocol="${strandedness}_${librarytype}_${pairing}_${readlength}"
         echo -e "${sample_name}\t${protocol}\t${disease_code}\t${disease_name}\t${color}\t${projects}\tFresh/Frozen" | sed 's/"//g' >> ${covariates_file}
      done

      # Get metadata for DFCI files
      dfci_metadata=$(echo ${dfci_ids} | xargs python3 /stjude/bin/bulk_describe.py -p project-F5444K89PZxXjBqVJ3Pp79B4 --ids)
      echo ${dfci_metadata} > dfci_metadata.json
   fi

   # Combine reference data with DFCI metadata
   # The "-s" argument to jq takes the input files and puts the objects in an array.
   # so .[0] + .[1] concatenates the first and second elements together, in this case,
   # both elements are arrays. So the output is a concatentated array.
   if [ "${include_pdx}" == "true" ]
   then
      jq -s add metadata.json dfci_metadata.json  > combined_metadata.json
   else
      cp metadata.json combined_metadata.json
   fi


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
         id=${input_counts[$i]}
         id=$(echo $id | jq -r '.["$dnanexus_link"]')

         j=$(dx describe --json $id)
         file_name=$(echo $j | jq -r '.name')
         sample_name=$(get_sample_name "$j")
         if [[ $sample_name == "null" ]]
         then
            sample_name=$(echo ${file_name} | cut -f 1 -d'.' )
         fi

         strandedness=$(get_property "$j" 'strandedness')
         if [[ $strandedness == "null" ]] && [[ ! -z $all_strandedness ]]
         then
            strandedness="${all_strandedness}"
         fi

         librarytype=$(get_property "$j" 'library_type')
         if [[ $librarytype == "null" ]] && [[ ! -z $all_library_type ]]
         then
            librarytype="${all_library_type}"
         fi

         readlength=$(get_property "$j" 'read_length')
         if [[ $readlength == "null" ]] && [[ ! -z $all_read_length ]]
         then
            readlength="${all_read_length}"
         fi

         pairing=$(get_property "$j" 'pairing')
         if [[ $pairing == "null" ]] && [[ ! -z $all_pairing ]]
         then
            pairing="${all_pairing}"
         fi

         disease_code=$sample_name
         protocol="${strandedness}_${librarytype}_${pairing}_${readlength}"

         # Check for missing properties on input file
         if [[ $protocol == "null_null_null_null" ]]
         then
            echo "{\"error\": {\"type\": \"AppError\", \"message\": \"Input file missing properties: $file_name\"}}" > job_error.json
            exit 1
         fi
         if [[ $strandedness == "null" ]]
         then
            echo "{\"error\": {\"type\": \"AppError\", \"message\": \"Input file missing property (strandedness): $file_name\"}}" > job_error.json
            exit 1
         fi
         if [[ $librarytype == "null" ]]
         then
            echo "{\"error\": {\"type\": \"AppError\", \"message\": \"Input file missing property (library_type): $file_name\"}}" > job_error.json
            exit 1
         fi
         if [[ $pairing == "null" ]]
         then
            echo "{\"error\": {\"type\": \"AppError\", \"message\": \"Input file missing property (pairing): $file_name\"}}" > job_error.json
            exit 1
         fi
         if [[ $readlength == "null" ]]
         then
            echo "{\"error\": {\"type\": \"AppError\", \"message\": \"Input file missing property (read_length): $file_name\"}}" > job_error.json
            exit 1
         fi

         echo -e "${sample_name}\t${protocol}\t${disease_code}\t\t\tinput\t" | sed 's/"//g' >> ${covariates_file}
         in_arg="$in_arg --input-sample $sample_name"
         echo "Adding input sample: $sample_name, $protocol, $disease_code"
      done
      echo ${in_arg}
   fi

   # Check covariates for single sample batch
   covariate_counts=$(cut -f 2 ${covariates_file} | sort |grep -v "Protocol" | awk '{a[$1]++}END{for (i in a) print i,a[i] | "sort"}' OFS="\t")
   echo "Covariates:"
   echo "$covariate_counts"
   smallest=$(echo "${covariate_counts}" | sort -k 2,2n | head -n 1)
   echo "smallest: ${smallest}"
   if [[ $(echo "${smallest}" | cut -f 2  )  == "1" ]]
   then
         co=$(echo "${smallest}" | cut -f 1)
         echo "{\"error\": {\"type\": \"AppError\", \"message\": \"A covariate batch has a single sample. This is unsupported for batch correction.: ${co}\"}}" > job_error.json
         exit 1
   fi

   # Fetch gene excludelist
   echo ""
   echo "  [*] Downloading gene excludelist ..."
   dx download -o $local_reference_dir/gene.excludelist.tsv project-F5444K89PZxXjBqVJ3Pp79B4:file-Fk84jFj97xxp1jxP9Zp6JJF4

   # Fetch Gencode
   echo ""
   echo "  [*] Downloading gencode ..."
   dx download -o $local_reference_dir/gencode.v31.annotation.gtf.gz -r project-F5444K89PZxXjBqVJ3Pp79B4:/pipeline/M2A/gencode.v31.annotation.gtf.gz

   # Run interactive t-SNE
   echo ""
   echo "  [*] Running t-SNE ..."
   tissue_arg=
   if [[ "$tumor_type" != "All" ]]
   then
      tissue_arg="--tissue-type \"${tumor_type}\""
   fi

   # Fetch Gene List
   gene_list_arg=
   if [[ ! -z $gene_list ]]
   then
      dx download -o $local_reference_dir/gene_list.txt "${gene_list}"
      gene_list_arg="--gene-list ${container_reference_dir}/gene_list.txt"
   fi

   echo "docker run -v $local_data_dir:$container_data_dir -v $local_reference_dir:$container_reference_dir -v $local_output_dir:$container_output_dir ghcr.io/stjudecloud/expression-classification:EXPRESSION_VERSION bash -c \"cd $container_output_dir && itsne-main --debug-rscript -b $container_reference_dir/gene.excludelist.tsv -g $container_reference_dir/gencode.v31.annotation.gtf.gz -c $container_data_dir/covariates.txt -o $container_output_dir/${output_name} ${in_arg} ${infile_arg} $container_data_dir/reference_counts/*.txt --save-data ${tissue_arg} ${gene_list_arg}\""

   docker run -v $local_data_dir:$container_data_dir -v $local_reference_dir:$container_reference_dir -v $local_output_dir:$container_output_dir ghcr.io/stjudecloud/expression-classification:EXPRESSION_VERSION bash -c "cd $container_output_dir && itsne-main --debug-rscript -b $container_reference_dir/gene.excludelist.tsv -g $container_reference_dir/gencode.v31.annotation.gtf.gz -c $container_data_dir/covariates.txt -o $container_output_dir/${output_name} ${in_arg} ${infile_arg} $container_data_dir/reference_counts/*.txt --save-data ${tissue_arg} ${gene_list_arg}"

   cp combined_metadata.json $local_output_dir
   cp /stjude/metadata/Subtype_Groupings_for_tSNE.csv $local_output_dir
   docker run -v $local_output_dir:$container_output_dir -v /stjude/bin:/stjude/bin ghcr.io/stjudecloud/expression-classification:EXPRESSION_VERSION bash -c "cd $container_output_dir && python /stjude/bin/generate_plot.py --tsne-file $container_output_dir/tsne.txt --metadata-file $container_output_dir/combined_metadata.json --subtype-file $container_output_dir/Subtype_Groupings_for_tSNE.csv $tissue_arg"

   if [ ${#input_counts[@]} -gt 0 ]
   then
      docker run -v $local_output_dir:$container_output_dir -v /stjude/bin:/stjude/bin ghcr.io/stjudecloud/expression-classification:EXPRESSION_VERSION bash -c "cd $container_output_dir && python /stjude/bin/neighbors.py tsne.txt neighbors.tsv"
      mv $local_output_dir/neighbors.tsv $local_output_dir/${output_prefix}.neighbors.tsv
      neighbors_file=$(dx upload $local_output_dir/${output_prefix}.neighbors.tsv --brief)
      dx-jobutil-add-output neighbors "$neighbors_file" --class=file
   fi

   # Upload output
   mv $local_output_dir/tsne_pp.html $local_output_dir/${output_name}
   tsne_plot=$(dx upload $local_output_dir/${output_name} --brief)
   dx-jobutil-add-output tsne_plot "$tsne_plot" --class=file
   mv $local_output_dir/tsne.txt $local_output_dir/${output_prefix}.tsne.txt
   tsne_matrix=$(dx upload $local_output_dir/${output_prefix}.tsne.txt --brief)
   dx-jobutil-add-output tsne_matrix "$tsne_matrix" --class=file
   dx-jobutil-add-output trustworthiness_score "$(cat $local_output_dir/trustworthiness.txt)" --class=string
   if [ -e $local_output_dir/gene_list.txt ]
   then
      mv $local_output_dir/gene_list.txt $local_output_dir/${output_prefix}.gene_list.txt
      gene_list=$(dx upload $local_output_dir/${output_prefix}.gene_list.txt --brief)
      dx-jobutil-add-output gene_list "$gene_list" --class=file
   fi
   if [ ${intermediate} ]
   then
      mv $local_output_dir/${intermediate}.txt $local_output_dir/${output_prefix}.${intermediate}.txt
      intermediate_file=$(dx upload $local_output_dir/${output_prefix}.${intermediate}.txt --brief)
      dx-jobutil-add-output intermediate_output "$intermediate_file" --class=file
   fi
}
