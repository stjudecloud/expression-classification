<!-- dx-header -->
# Interactive t-SNE (DNAnexus Platform App)

Produce a t-SNE plot that highlights input sample(s) against a cohort of existing samples.

This is the source code for an app that runs on the DNAnexus Platform.
For more information about how to run or modify it, see
https://wiki.dnanexus.com/.
<!-- /dx-header -->

<!-- Insert a description of your app here -->
This app assumes you have a reference dataset of counts files that are annotated with metadata about the samples. A reference dataset can be obtained by requesting it through [St. Jude Cloud](https://platform.stjude.cloud/data/publications?publication_accession=SJC-PB-1020).

---

## Parameters

Reference HTSeq count files should be provided using the `reference_counts` parameter. These can be obtained by requesting data through St. Jude Cloud Genomics Platform.

This app optionally takes input samples(s) and a tissue type. It then plots the t-Distributed Stochastic Neighbor Embedding (t-SNE) for the input samples relative to a cohort of samples.  

Input samples must have the following properties: "sample_name", "sj_diseases", "sample_type", "library_type", "read_length", "strandedness", and "pairing".  

* sample_name: should match the first part of the filename up to the first period (".")
* sj_diseases: a label for the category of the sample
* sample_type: diagnosis, relapse, etc.
* library_type: PolyA or Total
* read_length: read length in base pair
* strandedness: Stranded-Forward, Stranded-Reverse, or Unstranded
* pairing: Paired-end or Single-end

---

## WARNING

The RNA-Seq Expression Classification pipeline reference data is based on GRCh38 aligned, Gencode v31 annotated samples from fresh, frozen tissue samples. It has not been evaluated for samples that do not meet this criteria.

The RNA-Seq Expression Classification pipeline reference data uses sequencing data from fresh, frozen tissue samples. It has not been evaluated for use with sequencing data generated from formalin-fixed paraffin-embedded (FFPE) specimens.

If running the count-based RNA-Seq Expression Classification pipeline, alignment must be done against the GRCh38_no_alt reference. It should use parameters as specified in our RNA-seq workflow to minimize any discrepancies caused by differing alignment specification.

If running the count-based RNA-Seq Expression Classification pipeline, feature counts should be generated with htseq-count as described in our RNA-seq workflow. This pipeline uses Gencode v31 annotations.

Batch correction requires a minimum of two samples per batch to run properly. Introducing a single sample batch by adding an input sample with a unique protocol will cause unexpected results.

---

## Additional Information

Full documentation for running can be found in the [St. Jude Cloud docs](https://www.stjude.cloud/docs/guides/genomics-platform/analyzing-data/rnaseq-expression-classification/).

The t-SNE method is described in Laurens van der Maatens and Geoffrey Hinton. Visualizing Data using t-SNE. (2008)

<!--
TODO: This app directory was automatically generated by dx-app-wizard;
please edit this Readme.md file to include essential documentation about
your app that would be helpful to users. (Also see the
Readme.developer.md.) Once you're done, you can remove these TODO
comments.

For more info, see https://wiki.dnanexus.com/Developer-Portal.
-->
