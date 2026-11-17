# Public bulk RNA-seq: GSE271812

## Recovered source

`GSE271812_original_downstream_analysis_REVIEW.R` contains **3,565 lines of actual user-provided analysis commands** from the earlier `new.R` file, recovered from the chat file library. Another pasted version dated 2026-11-08 was checked and found text-identical after normalizing line endings. This is an **original downstream working transcript**, not a tested standalone pipeline or a newly invented reconstruction. The only functional edit was disabling a nonessential `.Rhistory` inspection snippet; a provenance header was added.

## Analyses found in the script

- Spearman correlation of `Ism1` and other genes in influenza-infected lung samples, with multiple explorations of day adjustment;
- GO biological process / KEGG enrichment with `gseGO` and related GSEA functions;
- comparisons and displays involving selected inflammation and innate immune pathways;
- GSVA pathway scores and Spearman plots comparing `Ism1` VST expression against pathway scores;
- day-4 subsets, DESeq2/VST preparation and related GSVA plots;
- exports of CSV, PNG and PDF results.

**Important:** Some sections are exploratory duplicates, and early titles in this historical script say “Figure 3” even though current thesis figure numbering differs. Do not use titles as proof that each plot is in the submitted thesis.

## Required prior R objects / missing upstream steps

The script **starts** with:

```r
stopifnot(exists("vst.matrix.iav271812"), exists("metadata.timecourse271812"))
```

and later refers to objects such as `count.matrix.primary271812`, `meta.gsva.271812`, and `gsva.plot.data.271812` whose original creation/loading code is not present in this recovered script. This means that **importing GSE271812 counts, constructing the complete design metadata, generating the initial VST expression matrix, and some GSVA calculations remain to be located or reconstructed and checked against the user's real data**. Do not describe this version as fully runnable from a fresh R session.

## Next checks before public release

1. Identify the earlier script, RData/RDS, or notebook that creates the prerequisite objects; use genuine data and sample metadata, do not invent analysis steps.
2. Check condition/timepoint definitions and gene IDs against the Methods and Figure 2 legends.
3. Verify the precise scoring/ranking inputs, correlation method, multiplicity handling, and GSEA cutoffs for the *final figures*, rather than copying settings from exploratory sections.
4. Test a clean re-run under recorded R and package versions.
5. Obtain supervisor approval before publishing any scripts associated with unpublished laboratory results.

The raw data and R objects are **not included** in this repository. No credentials were copied into this folder.
