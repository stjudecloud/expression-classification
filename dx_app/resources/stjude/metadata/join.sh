# Add "key" column as join of project and sample name to allow unique matches
# CLINPILOT -> Clinical Pilot in SAMPLE_INFO_August2019_Working_V99_03_August_FREEZE_V4_slim.csv
# Set 30 RTCG samples to EXCLUDE due to missing metadata: SJAML030612_D1 SJAML031051_D1,SJAML031106_D1,SJAML031120_D1,SJAML031121_D1,SJBALL031076_D1,SJBALL031079_D1,SJBALL031115_D1,SJBALL031116_D1,SJBALL031128_D1,SJBALL031144_D1,SJBT031135_D1,SJBT031136_D1,SJHGG030775_D2,SJHGG031097_D1,SJHGG031140_D1,SJHM030441_D2,SJHM031090_D1,SJMB031110_D1,SJNBL030452_D1,SJNBL031046_D2,SJNBL031145_D1,SJOS031125_D1,SJOS031130_D1,SJRHB031084_D1,SJRHB031117_D1,SJST031108_D1,SJST031143_D1,SJTALL031107_D1,SJWLM030577_D1,
# Updated APED -> AEPD in SAMPLE_INFO_August2019_Working_V99_03_August_FREEZE_V4_slim.csv
# Updated MSCERMS -> MSCERHB in SAMPLE_INFO_August2019_Working_V99_03_August_FREEZE_V4_slim.csv

csvjoin -c 1  --left PCGP_RTCG_G4K_ClinGen_RNAseqInfo_Aug_28_2019_FREEZE_V4_project.csv SAMPLE_INFO_August2019_Working_V99_03_August_FREEZE_V4_slim.csv > combined1.csv
csvjoin -c 30,2 --left combined1.csv SAMPLE_INFO_August2019_Working_V99_LOOKUP_V4_slim.csv > combined2.csv
csvcut -c SampleID,Sample_Preparation,Sequencing_Machine,Library_Selection,Sequencing_Approach,ReadLength_bp,"Predicted Strandness via InferExperiment",Dataset,"TumorCategory NEW","Classification_tSNE_tag (Manuscript)","Classification_tSNE_color (Manuscript)",KeepOrExclude combined2.csv |grep -v "EXCLUDE" > combined.csv
