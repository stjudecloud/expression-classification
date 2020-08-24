# Add "key" column as join of project and sample name to allow unique matches paper_vs_database_diagnosis_v4_normalized_AlexUpdatesV2 and PCGP_RTCG_G4K_ClinGen_RNAseqInfo_Aug_28_2019_FREEZE_V6
# Set 30 RTCG samples to EXCLUDE due to missing metadata: SJAML030612_D1 SJAML031051_D1,SJAML031106_D1,SJAML031120_D1,SJAML031121_D1,SJBALL031076_D1,SJBALL031079_D1,SJBALL031115_D1,SJBALL031116_D1,SJBALL031128_D1,SJBALL031144_D1,SJBT031135_D1,SJBT031136_D1,SJHGG030775_D2,SJHGG031097_D1,SJHGG031140_D1,SJHM030441_D2,SJHM031090_D1,SJMB031110_D1,SJNBL030452_D1,SJNBL031046_D2,SJNBL031145_D1,SJOS031125_D1,SJOS031130_D1,SJRHB031084_D1,SJRHB031117_D1,SJST031108_D1,SJST031143_D1,SJTALL031107_D1,SJWLM030577_D1,

csvjoin -c 7,3 paper_vs_database_diagnosis_v4_normalized_AlexUpdatesV2.csv paper_vs_database_diagnosis_v4_normalized_AlexUpdatesV2_LOOKUP_tSNE.csv > combined1.csv 
csvjoin -c 1,1 combined1.csv PCGP_RTCG_G4K_ClinGen_RNAseqInfo_Aug_28_2019_FREEZE_V6.csv > combined2.csv  
csvcut -c SampleID,Sample_Preparation,Sequencing_Machine,Library_Selection,Sequencing_Approach,ReadLength_bp,"Predicted Strandness via InferExperiment",Dataset,"TumorCategory NEW","St. Jude Diagnosis ID NEW","Classification_tSNE_color (Manuscript)",KeepOrExclude combined2.csv |grep -v "EXCLUDE" > combined.csv
