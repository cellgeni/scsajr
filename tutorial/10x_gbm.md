- [1 Prepare inputs](#prepare-inputs)
  - [1.1 Install packages (run once)](#install-packages-run-once)
  - [1.2 Load packages](#load-packages)
  - [1.3 Set working directory](#set-working-directory)
  - [1.4 Load reference](#load-reference)
  - [1.5 Load input data](#load-input-data)
- [2 Inspect pipeline output](#inspect-pipeline-output)
  - [2.1 Load Pseudobulk Splicing Data](#load-pseudobulk-splicing-data)
    - [2.1.1 Row: Information about segments
      (features)](#row-information-about-segments-features)
    - [2.1.2 information about pseudobulks (per
      sample/celltype)](#information-about-pseudobulks-per-samplecelltype)
    - [2.1.3 assays](#assays)
    - [2.1.4 metadata: postprocessing
      results](#metadata-postprocessing-results)
- [3 Find differentially spliced
  segments](#find-differentially-spliced-segments)
  - [3.1 Results comparing all celltypes simultaneously
    $`\text{fdr} < 0.05`$](#results-comparing-all-celltypes-simultaneously-textfdr-0.05)
  - [3.2 Results of comparing each celltype against the rest are stored
    as three (pv, fdr, dpsi) $`\text{segment} * \text{celltype}`$
    matrices](#results-of-comparing-each-celltype-against-the-rest-are-stored-as-three-pv-fdr-dpsi-textsegment-textcelltype-matrices)
  - [3.3 Marker segment selection](#marker-segment-selection)
  - [3.4 Load gene expression data](#load-gene-expression-data)
- [4 Convert Seurat to SummarizedExperiment & pseudobulk
  expression](#convert-seurat-to-summarizedexperiment-pseudobulk-expression)
  - [4.1 Visualise markers](#visualise-markers)
  - [4.2 Coverage plots](#coverage-plots)
  - [4.3 Differential splicing](#differential-splicing)
  - [4.4 domain enrichment analyses](#domain-enrichment-analyses)
  - [4.5 re-pseudobulk](#re-pseudobulk)

# 1 Prepare inputs

## 1.1 Install packages (run once)

``` r
# install.packages("BiocManager")
# BiocManager::install("EBImage")
# BiocManager::install("GenomicAlignments")
# BiocManager::install("clusterProfiler")
# devtools::install_github("cellgeni/visutils")
# devtools::install_github("mamarkevi/plotCoverage")
# devtools::install_github("cellgeni/scsajr")
```

## 1.2 Load packages

``` r
library(Seurat)
library(plotCoverage)
library(clusterProfiler)
library(Matrix)
library(SummarizedExperiment)
library(scsajr)
```

## 1.3 Set working directory

``` r
data_path <- "." # path to tutorial folder
setwd(data_path)
```

## 1.4 Load reference

``` r
# gene_descr[dataframe]: gene descriptions
gene_descr_path <- "./nf-scsajr/ref/human_2020A/functional_annotation/gene.descr.rds"
gene_descr <- readRDS(gene_descr_path)
str(gene_descr)
```

    ## 'data.frame':    36601 obs. of  8 variables:
    ##  $ ens_id: chr  "ENSG00000243485" "ENSG00000237613" "ENSG00000186092" "ENSG00000238009" ...
    ##  $ descr : chr  "MIR1302-2 host gene" "family with sequence similarity 138 member A" "olfactory receptor family 4 subfamily F member 5" "novel transcript" ...
    ##  $ chr   : chr  "1" "1" "1" "1" ...
    ##  $ start : int  29554 34554 65419 89295 89551 139790 141474 160446 266855 358857 ...
    ##  $ end   : int  31109 36081 71585 133723 91105 140339 173862 161525 268655 366052 ...
    ##  $ strand: int  1 -1 1 -1 -1 -1 -1 1 1 1 ...
    ##  $ name  : chr  "MIR1302-2HG" "FAM138A" "OR4F5" "AL627309.1" ...
    ##  $ type  : chr  "lncRNA" "lncRNA" "protein_coding" "lncRNA" ...

``` r
# View(gene_descr)
```

``` r
# gtf[dataframe]: GTF file with gene annotations
gtf_path <- "./nf-scsajr/ref/human_2020A/gtf.rds"
gtf <- readRDS(gtf_path)
str(gtf)
```

    ## 'data.frame':    2765969 obs. of  26 variables:
    ##  $ chr_id                  : chr  "1" "1" "1" "1" ...
    ##  $ type                    : chr  "HAVANA" "HAVANA" "HAVANA" "HAVANA" ...
    ##  $ feature                 : chr  "gene" "transcript" "exon" "exon" ...
    ##  $ start                   : int  29554 29554 29554 30564 30976 30267 30267 30976 34554 34554 ...
    ##  $ stop                    : int  31109 31097 30039 30667 31097 31109 30667 31109 36081 36081 ...
    ##  $ strand                  : chr  "+" "+" "+" "+" ...
    ##  $ gene_id                 : chr  "ENSG00000243485" "ENSG00000243485" "ENSG00000243485" "ENSG00000243485" ...
    ##  $ gene_version            : chr  "5" "5" "5" "5" ...
    ##  $ gene_type               : chr  "lncRNA" "lncRNA" "lncRNA" "lncRNA" ...
    ##  $ gene_name               : chr  "MIR1302-2HG" "MIR1302-2HG" "MIR1302-2HG" "MIR1302-2HG" ...
    ##  $ level                   : chr  "2" "2" "2" "2" ...
    ##  $ hgnc_id                 : chr  "HGNC:52482" "HGNC:52482" "HGNC:52482" "HGNC:52482" ...
    ##  $ tag                     : chr  "ncRNA_host" "not_best_in_genome_evidence" "not_best_in_genome_evidence" "not_best_in_genome_evidence" ...
    ##  $ havana_gene             : chr  "OTTHUMG00000000959.2" "OTTHUMG00000000959.2" "OTTHUMG00000000959.2" "OTTHUMG00000000959.2" ...
    ##  $ transcript_id           : chr  NA "ENST00000473358" "ENST00000473358" "ENST00000473358" ...
    ##  $ transcript_version      : chr  NA "1" "1" "1" ...
    ##  $ transcript_type         : chr  NA "lncRNA" "lncRNA" "lncRNA" ...
    ##  $ transcript_name         : chr  NA "MIR1302-2HG-202" "MIR1302-2HG-202" "MIR1302-2HG-202" ...
    ##  $ transcript_support_level: chr  NA "5" "5" "5" ...
    ##  $ havana_transcript       : chr  NA "OTTHUMT00000002840.1" "OTTHUMT00000002840.1" "OTTHUMT00000002840.1" ...
    ##  $ exon_number             : chr  NA NA "1" "2" ...
    ##  $ exon_id                 : chr  NA NA "ENSE00001947070" "ENSE00001922571" ...
    ##  $ exon_version            : chr  NA NA "1" "1" ...
    ##  $ protein_id              : chr  NA NA NA NA ...
    ##  $ ccdsid                  : chr  NA NA NA NA ...
    ##  $ ont                     : chr  NA NA NA NA ...

``` r
# View(gtf)  # Too large to view
```

## 1.5 Load input data

``` r
# Sample -> BAM file (scRNA data mapped to reference genome)
# echo -e "sample\t$(pwd)/alignment.bam" > samples.tsv
samples <- read.table("./input/samples.tsv", col.names = c("sample_id", "bam_path"))

# Sample -> barcode, celltype
# awk -F, 'BEGIN { OFS="\t" }; NR>1 {print  "sample",$1,$2}' analysis/clustering/graphclust/clusters.csv > barcodes.tsv
barcodes <- read.table("./input/barcodes.tsv", col.names = c("sample_id", "barcode", "celltype"))
rownames(barcodes) <- paste0(barcodes$sample_id, "|", barcodes$barcode)

# Definition of splicing segments
segs_path <- "./nf-scsajr/ref/human_2020A/segments.csv"
segs <- read.csv(segs_path, row.names = 1)
```

------------------------------------------------------------------------

# 2 Inspect pipeline output

## 2.1 Load Pseudobulk Splicing Data

``` r
# pbas.rds: Pseudobulk AS segments data for all segments
# pb_as_filtered.rds [SummerizedExperiments]: Pseudobulk AS segments data filtered by coverage
pbasf_path <- "./nf-scsajr_output/rds/pb_as_filtered.rds"
pbasf <- readRDS(pbasf_path)

slotNames(pbasf)
```

    ## [1] "rowRanges"       "colData"         "assays"          "NAMES"          
    ## [5] "elementMetadata" "metadata"

| Column | Meaning |
|----|----|
| `Rows` | Splicing segments |
| `Cols` | (sample, cell-type) pseudobulks |
| `Assays` | `i`: inclusion junction, `e`: exclusion junction, `psi`: percent spliced |
| `Metadata` | Differential-splicing test results |

<br><br>

### 2.1.1 Row: Information about segments (features)

| Column | Meaning |
|----|----|
| `feature_id` | unique ID from Ensembl gene ID “ENSG00000188976.s7” |
| `gene_id` | Ensembl gene ID (e.g. “ENSG00000188976”) |
| `type` | segment type (e.g. “EXN” for exon, “ALT” for alternative region) |
| `position` | location within transcript (e.g. “INTERNAL”) |
| `sites` | splice‐site pattern (e.g. “ad” = acceptor–donor) |
| `length` | width of the segment in bp |
| `is_exon` | logical: is this a cassette exon? |
| `cod`, `cod.gene` | coding status of the segment (within CDS) |
| `ncell` | number of single cells with coverage ≥10 for this segment |
| `nna` | how many pseudobulks had `NA` for psi (insufficient reads) |
| `sd` | standard deviation of ψ across pseudobulks |

``` r
SummarizedExperiment::rowRanges(pbasf)[1:3]
```

    ## GRanges object with 3 ranges and 12 metadata columns:
    ##                      seqnames        ranges strand |         feature_id
    ##                         <Rle>     <IRanges>  <Rle> |        <character>
    ##   ENSG00000188976.s7     chr1 946402-946545      - | ENSG00000188976.s7
    ##   ENSG00000188976.s9     chr1 948131-948232      - | ENSG00000188976.s9
    ##   ENSG00000188290.s4     chr1 999692-999787      - | ENSG00000188290.s4
    ##                              gene_id        type    position       sites
    ##                          <character> <character> <character> <character>
    ##   ENSG00000188976.s7 ENSG00000188976         EXN    INTERNAL          ad
    ##   ENSG00000188976.s9 ENSG00000188976         EXN    INTERNAL          ad
    ##   ENSG00000188290.s4 ENSG00000188290         ALT    INTERNAL          ad
    ##                         length   is_exon         cod  cod.gene     ncell
    ##                      <integer> <logical> <character> <logical> <numeric>
    ##   ENSG00000188976.s7       144      TRUE           c      TRUE         0
    ##   ENSG00000188976.s9       102      TRUE           c      TRUE         0
    ##   ENSG00000188290.s4        96      TRUE           c      TRUE       282
    ##                            nna        sd
    ##                      <integer> <numeric>
    ##   ENSG00000188976.s7         5  0.108319
    ##   ENSG00000188976.s9         5  0.120148
    ##   ENSG00000188290.s4         8  0.172262
    ##   -------
    ##   seqinfo: 40 sequences from an unspecified genome; no seqlengths

``` r
# View(pbasf@rowRanges)
```

<br><br>

### 2.1.2 information about pseudobulks (per sample/celltype)

| Column      | Meaning                              |
|-------------|--------------------------------------|
| `sample_id` | original sample ID                   |
| `celltype`  | the celltype (cluster)               |
| `ncells`    | cells aggregated into the pseudobulk |
| `strand`    | strand orientation of the coverage   |

``` r
pbasf@colData
```

    ## DataFrame with 10 rows and 4 columns
    ##             sample_id    celltype    ncells    strand
    ##           <character> <character> <numeric> <integer>
    ## sample$1       sample           1       716        -1
    ## sample$10      sample          10       194        -1
    ## sample$2       sample           2       550        -1
    ## sample$3       sample           3       538        -1
    ## sample$4       sample           4       459        -1
    ## sample$5       sample           5       369        -1
    ## sample$6       sample           6       350        -1
    ## sample$7       sample           7       341        -1
    ## sample$8       sample           8       301        -1
    ## sample$9       sample           9       285        -1

<br><br>

### 2.1.3 assays

| Assay name | Contents                                                 |
|------------|----------------------------------------------------------|
| `i`        | inclusion junction counts                                |
| `e`        | exclusion junction counts                                |
| `psi`      | percent-spliced-in $`\frac{i}{(i+e)}`$; NA if $`i+e<10`$ |

``` r
pbasf@assays
```

    ## An object of class "SimpleAssays"
    ## Slot "data":
    ## List of length 3
    ## names(3): i e psi

``` r
SummarizedExperiment::assay(pbasf, "i")[1:5, 1:4]
```

    ##                    sample$1 sample$10 sample$2 sample$3
    ## ENSG00000188976.s7        9         0       15       13
    ## ENSG00000188976.s9       15         0       14        6
    ## ENSG00000188290.s4      940        62     1263     1478
    ## ENSG00000186891.s4       14         2       13        3
    ## ENSG00000186891.s5       23         0       18       14

``` r
SummarizedExperiment::assay(pbasf, "e")[1:5, 1:4]
```

    ##                    sample$1 sample$10 sample$2 sample$3
    ## ENSG00000188976.s7        0         0        0        0
    ## ENSG00000188976.s9        0         0        0        0
    ## ENSG00000188290.s4     1444       116      336      626
    ## ENSG00000186891.s4        5         0        0        0
    ## ENSG00000186891.s5       11         2        8        1

``` r
SummarizedExperiment::assay(pbasf, "psi")[1:5, 1:4]
```

    ##                     sample$1 sample$10  sample$2  sample$3
    ## ENSG00000188976.s7        NA        NA 1.0000000 1.0000000
    ## ENSG00000188976.s9 1.0000000        NA 1.0000000        NA
    ## ENSG00000188290.s4 0.3942953 0.3483146 0.7898687 0.7024715
    ## ENSG00000186891.s4 0.7368421        NA 1.0000000        NA
    ## ENSG00000186891.s5 0.6764706        NA 0.6923077 0.9333333

<br><br>

### 2.1.4 metadata: postprocessing results

- `all_celltype_test`
  - `group`/`group_fdr` — test statistics & FDR
  - `low_state`, `high_state` — celltype with lowest and highest
  - `dpsi` — difference in psi between groups
- `markers`
  - `pv` — p-value per segment and celltype
  - `fdr` — FDR per segment and celltype
  - `dpsi` — difference in psi per segment and celltype
- `go` — GO enrichment results
- `ipro` — domain enrichment results

``` r
names(pbasf@metadata)
```

    ## [1] "all_celltype_test" "markers"           "go"               
    ## [4] "ipro"

``` r
head(S4Vectors::metadata(pbasf)$all_celltype_test)
```

    ##                    overdispersion         group     group_fdr low_state
    ## ENSG00000188976.s7             NA  1.734812e-01  1.804119e-01      <NA>
    ## ENSG00000188976.s9             NA  6.041153e-02  6.969426e-02      <NA>
    ## ENSG00000188290.s4             NA 8.452251e-195 4.101757e-192        10
    ## ENSG00000186891.s4             NA  6.930312e-02  7.860525e-02      <NA>
    ## ENSG00000186891.s5             NA  3.593649e-05  1.227514e-04      <NA>
    ## ENSG00000186891.s6             NA  1.912381e-05  7.084362e-05         1
    ##                    high_state       dpsi
    ## ENSG00000188976.s7       <NA>         NA
    ## ENSG00000188976.s9       <NA>         NA
    ## ENSG00000188290.s4          2 0.44155406
    ## ENSG00000186891.s4       <NA>         NA
    ## ENSG00000186891.s5       <NA>         NA
    ## ENSG00000186891.s6          5 0.06025563

``` r
str(S4Vectors::metadata(pbasf)$markers)
```

    ## List of 3
    ##  $ pv  : num [1:6794, 1:10] 4.87e-01 2.76e-01 2.95e-04 6.73e-05 8.24e-01 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:6794] "ENSG00000188976.s7" "ENSG00000188976.s9" "ENSG00000188290.s4" "ENSG00000186891.s4" ...
    ##   .. ..$ : chr [1:10] "1" "10" "2" "3" ...
    ##  $ fdr : num [1:6794, 1:10] 0.84899 0.73889 0.01398 0.00478 0.93365 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:6794] "ENSG00000188976.s7" "ENSG00000188976.s9" "ENSG00000188290.s4" "ENSG00000186891.s4" ...
    ##   .. ..$ : chr [1:10] "1" "10" "2" "3" ...
    ##  $ dpsi: num [1:6794, 1:10] NaN 0.0976 -0.2174 -0.2632 -0.1473 ...
    ##   ..- attr(*, "dimnames")=List of 2
    ##   .. ..$ : chr [1:6794] "ENSG00000188976.s7" "ENSG00000188976.s9" "ENSG00000188290.s4" "ENSG00000186891.s4" ...
    ##   .. ..$ : chr [1:10] "1" "10" "2" "3" ...

``` r
str(S4Vectors::metadata(pbasf)$go)
```

    ## Formal class 'compareClusterResult' [package "DOSE"] with 10 slots
    ##   ..@ compareClusterResult:'data.frame': 33 obs. of  14 variables:
    ##   .. ..$ Cluster       : Factor w/ 11 levels "1","10","2","3",..: 2 2 2 2 2 5 5 6 6 6 ...
    ##   .. ..$ ONTOLOGY      : chr [1:33] "BP" "BP" "BP" "BP" ...
    ##   .. ..$ ID            : chr [1:33] "GO:0032755" "GO:0032635" "GO:0032675" "GO:0072567" ...
    ##   .. ..$ Description   : chr [1:33] "positive regulation of interleukin-6 production" "interleukin-6 production" "regulation of interleukin-6 production" "chemokine (C-X-C motif) ligand 2 production" ...
    ##   .. ..$ GeneRatio     : chr [1:33] "3/11" "3/11" "3/11" "2/11" ...
    ##   .. ..$ BgRatio       : chr [1:33] "56/8527" "83/8527" "83/8527" "15/8527" ...
    ##   .. ..$ RichFactor    : num [1:33] 0.0536 0.0361 0.0361 0.1333 0.1333 ...
    ##   .. ..$ FoldEnrichment: num [1:33] 41.5 28 28 103.4 103.4 ...
    ##   .. ..$ zScore        : num [1:33] 10.94 8.89 8.89 14.26 14.26 ...
    ##   .. ..$ pvalue        : num [1:33] 4.27e-05 1.39e-04 1.39e-04 1.57e-04 1.57e-04 ...
    ##   .. ..$ p.adjust      : num [1:33] 0.0337 0.0337 0.0337 0.0337 0.0337 ...
    ##   .. ..$ qvalue        : num [1:33] 0.0238 0.0238 0.0238 0.0238 0.0238 ...
    ##   .. ..$ geneID        : chr [1:33] "ENSG00000276600/ENSG00000197971/ENSG00000105329" "ENSG00000276600/ENSG00000197971/ENSG00000105329" "ENSG00000276600/ENSG00000197971/ENSG00000105329" "ENSG00000197971/ENSG00000105329" ...
    ##   .. ..$ Count         : int [1:33] 3 3 3 2 2 11 12 5 4 3 ...
    ##   ..@ geneClusters        :List of 11
    ##   .. ..$ 1  : chr [1:184] "ENSG00000188290" "ENSG00000186891" "ENSG00000127483" "ENSG00000196182" ...
    ##   .. ..$ 10 : chr [1:11] "ENSG00000122218" "ENSG00000276600" "ENSG00000065600" "ENSG00000138381" ...
    ##   .. ..$ 2  : chr [1:61] "ENSG00000143147" "ENSG00000171163" "ENSG00000172456" "ENSG00000162692" ...
    ##   .. ..$ 3  : chr [1:90] "ENSG00000007923" "ENSG00000117500" "ENSG00000164011" "ENSG00000116212" ...
    ##   .. ..$ 4  : chr [1:161] "ENSG00000188529" "ENSG00000142733" "ENSG00000126106" "ENSG00000132128" ...
    ##   .. ..$ 5  : chr [1:36] "ENSG00000157873" "ENSG00000116747" "ENSG00000179818" "ENSG00000144118" ...
    ##   .. ..$ 6  : chr [1:116] "ENSG00000134668" "ENSG00000116885" "ENSG00000196517" "ENSG00000009307" ...
    ##   .. ..$ 7  : chr [1:6] "ENSG00000075945" "ENSG00000055917" "ENSG00000176623" "ENSG00000149541" ...
    ##   .. ..$ 8  : chr [1:157] "ENSG00000178922" "ENSG00000134222" "ENSG00000065600" "ENSG00000143748" ...
    ##   .. ..$ 9  : chr [1:191] "ENSG00000221978" "ENSG00000236963" "ENSG00000116350" "ENSG00000066185" ...
    ##   .. ..$ all: chr [1:466] "ENSG00000188290" "ENSG00000235098" "ENSG00000116670" "ENSG00000077549" ...
    ##   ..@ fun                 : chr "enrichGO"
    ##   ..@ gene2Symbol         : chr(0) 
    ##   ..@ keytype             : chr "ENSEMBL"
    ##   ..@ readable            : logi FALSE
    ##   ..@ .call               : language compareCluster(geneClusters = gids, fun = "enrichGO", universe = gene_uni,      pAdjustMethod = "BH", ont = "ALL"| __truncated__
    ##   ..@ termsim             : num[0 , 0 ] 
    ##   ..@ method              : chr(0) 
    ##   ..@ dr                  : list()

``` r
str(S4Vectors::metadata(pbasf)$ipro)
```

    ## Formal class 'compareClusterResult' [package "DOSE"] with 10 slots
    ##   ..@ compareClusterResult:'data.frame': 297 obs. of  13 variables:
    ##   .. ..$ Cluster       : Factor w/ 11 levels "1","10","2","3",..: 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. ..$ ID            : chr [1:297] "IPR008532" "PF05670" "PF05833" "IPR031887" ...
    ##   .. ..$ Description   : chr [1:297] "NFACT, RNA-binding domain" "NFACT protein RNA binding domain" "NFACT N-terminal and middle domains" "Serologically defined colon cancer antigen 8" ...
    ##   .. ..$ GeneRatio     : chr [1:297] "2/92" "2/92" "2/92" "2/92" ...
    ##   .. ..$ BgRatio       : chr [1:297] "11/79860" "11/79860" "15/79860" "17/79860" ...
    ##   .. ..$ RichFactor    : num [1:297] 0.182 0.182 0.133 0.118 0.118 ...
    ##   .. ..$ FoldEnrichment: num [1:297] 158 158 116 102 102 ...
    ##   .. ..$ zScore        : num [1:297] 17.7 17.7 15.1 14.2 14.2 ...
    ##   .. ..$ pvalue        : num [1:297] 7.17e-05 7.17e-05 1.36e-04 1.77e-04 1.77e-04 ...
    ##   .. ..$ p.adjust      : num [1:297] 0.00494 0.00494 0.00494 0.00494 0.00494 ...
    ##   .. ..$ qvalue        : num [1:297] 0.0029 0.0029 0.0029 0.0029 0.0029 ...
    ##   .. ..$ geneID        : chr [1:297] "ENSG00000165525.s23/ENSG00000165525.s26" "ENSG00000165525.s23/ENSG00000165525.s26" "ENSG00000165525.s33/ENSG00000165525.s35" "ENSG00000054282.s9/ENSG00000054282.s19" ...
    ##   .. ..$ Count         : int [1:297] 2 2 2 2 2 2 2 2 2 2 ...
    ##   ..@ geneClusters        :List of 11
    ##   .. ..$ 1  : chr [1:155] "ENSG00000188290.s4" "ENSG00000186891.s4" "ENSG00000127483.s19" "ENSG00000196182.s17" ...
    ##   .. ..$ 10 : chr [1:10] "ENSG00000122218.s41" "ENSG00000276600.s7" "ENSG00000065600.s7" "ENSG00000138381.s4" ...
    ##   .. ..$ 2  : chr [1:56] "ENSG00000143147.s9" "ENSG00000171163.s7" "ENSG00000172456.s16" "ENSG00000153363.s3" ...
    ##   .. ..$ 3  : chr [1:82] "ENSG00000007923.s4" "ENSG00000117500.s2" "ENSG00000164011.s4" "ENSG00000116212.s12" ...
    ##   .. ..$ 4  : chr [1:143] "ENSG00000188529.s5" "ENSG00000142733.s7" "ENSG00000132128.s10" "ENSG00000121310.s3" ...
    ##   .. ..$ 5  : chr [1:33] "ENSG00000179818.s74" "ENSG00000144118.s7" "ENSG00000181722.s22" "ENSG00000138758.s7" ...
    ##   .. ..$ 6  : chr [1:102] "ENSG00000116885.s19" "ENSG00000196517.s20" "ENSG00000009307.s24" "ENSG00000204138.s5" ...
    ##   .. ..$ 7  : chr [1:5] "ENSG00000075945.s21" "ENSG00000055917.s21" "ENSG00000176623.s20" "ENSG00000149541.s8" ...
    ##   .. ..$ 8  : chr [1:147] "ENSG00000178922.s5" "ENSG00000134222.s2" "ENSG00000065600.s5" "ENSG00000143748.s32" ...
    ##   .. ..$ 9  : chr [1:183] "ENSG00000221978.s5" "ENSG00000236963.s6" "ENSG00000116350.s2" "ENSG00000066185.s9" ...
    ##   .. ..$ all: chr [1:415] "ENSG00000188290.s4" "ENSG00000235098.s2" "ENSG00000116670.s15" "ENSG00000077549.s2" ...
    ##   ..@ fun                 : chr "enricher"
    ##   ..@ gene2Symbol         : chr(0) 
    ##   ..@ keytype             : chr "UNKNOWN"
    ##   ..@ readable            : logi FALSE
    ##   ..@ .call               : language compareCluster(geneClusters = sids, fun = "enricher", universe = seg_uni,      pAdjustMethod = "BH", TERM2GENE = | __truncated__ ...
    ##   ..@ termsim             : num[0 , 0 ] 
    ##   ..@ method              : chr(0) 
    ##   ..@ dr                  : list()

<br><br>

# 3 Find differentially spliced segments

## 3.1 Results comparing all celltypes simultaneously $`\text{fdr} < 0.05`$

<!-- Need to understand statistics -->

``` r
subset(
  pbasf@metadata$all_celltype_test,
  !is.na(dpsi) & group_fdr < 0.05
)
```

## 3.2 Results of comparing each celltype against the rest are stored as three (pv, fdr, dpsi) $`\text{segment} * \text{celltype}`$ matrices

``` r
pbasf@metadata$markers$fdr[1:2, ]
```

    ##                            1 10         2         3         4         5
    ## ENSG00000188976.s7 0.8489872 NA 0.8183147 0.8279779 0.9691289 0.9540755
    ## ENSG00000188976.s9 0.7388893 NA 0.8063200 0.8777340 0.9179100        NA
    ##                            6         7         8         9
    ## ENSG00000188976.s7 0.8448498 0.9618651 0.1258921 0.9176569
    ## ENSG00000188976.s9 0.8204081        NA 0.3050156 0.9482446

``` r
pbasf@metadata$markers$dpsi[1:2, ]
```

    ##                             1  10          2          3            4   5
    ## ENSG00000188976.s7        NaN NaN 0.07601351 0.07601351  0.008445946 NaN
    ## ENSG00000188976.s9 0.09759358 NaN 0.09759358        NaN -0.049465241 NaN
    ##                             6   7          8   9
    ## ENSG00000188976.s7 0.07601351 NaN -0.2364865 NaN
    ## ENSG00000188976.s9 0.09759358 NaN -0.2433155 NaN

## 3.3 Marker segment selection

markers: $`\text{fdr} < 0.05`$, $`|\text{dpsi}| > 0.5`$, at most 2
segments per celltype

``` r
# select_markers only uses output of per-celltype tests
# let's find at most two marker segments per celltype, ensuring that each segment is reported only once
scsajr::select_markers(pbasf@metadata$markers, n = 2, dpsi_thr = 0.5, clean_duplicates = TRUE)
```

``` r
# select_all_markers also takes results of "all celltypes together" test
# it gives priority to markers (they are marked by is_marker field)

# select_all_markers relies on select_markers and select_markers_from_all_celltype_test
scsajr::select_all_markers(pbasf@metadata$markers, pbasf@metadata$all_celltype_test, dpsi_thr = 0.5, n = Inf)
```

## 3.4 Load gene expression data

This dataset is targeted, meaning that some neuronal genes were
enriched, but other genes are still present. The AS analysis was run on
all genes.

``` r
ge <- Seurat::Read10X_h5("./input/raw_feature_bc_matrix.h5", use.names = FALSE) # Load 10x h5 file
ge <- Seurat::CreateSeuratObject(ge) # Wrapping into Seurat object
```

``` r
ge <- ge[, barcodes$barcode] # Keep only cells annotated in barcodes.tsv
ge$barcode <- colnames(ge) # Add barcode as a column
colnames(ge) <- paste0("sample|", colnames(ge)) # Rename column names to include sample_id in barcodes
ge$celltype <- barcodes[colnames(ge), "celltype"] # Add pre-computed celltype labels from barcodes.tsv
```

``` r
# Just use standard Seurat pipeline to re-cluster cells
ge <- Seurat::NormalizeData(ge) # Log-normalise UMI counts
```

    ## Normalizing layer: counts

``` r
ge <- Seurat::FindVariableFeatures(ge) # Find highly variable features
```

    ## Finding variable features for layer counts

``` r
# VariableFeatures(ge)
ge <- Seurat::ScaleData(ge)
```

    ## Centering and scaling data matrix

``` r
ge <- Seurat::RunPCA(ge)
ge <- Seurat::FindNeighbors(ge)
ge <- Seurat::FindClusters(ge)
```

    ## Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    ## 
    ## Number of nodes: 4103
    ## Number of edges: 121867
    ## 
    ## Running Louvain algorithm...
    ## Maximum modularity in 10 random starts: 0.8781
    ## Number of communities: 17
    ## Elapsed time: 0 seconds

``` r
ge <- Seurat::RunUMAP(ge, dims = 1:30)
```

    ## Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    ## To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    ## This message will be shown once per session

``` r
Seurat::DimPlot(ge, reduction = "umap", group.by = c("celltype", "seurat_clusters"))
```

![](10x_gbm_files/figure-gfm/PCA%20and%20UMAP-1.png)<!-- -->

``` r
# Here we have just one sample, lets add sample_id to make code compatible with multiple sample design
ge$sample_id <- "sample"
```

# 4 Convert Seurat to SummarizedExperiment & pseudobulk expression

``` r
genes <- gene_descr[rownames(ge), ]
rowRanges <- GenomicRanges::GRanges(genes$chr,
  IRanges::IRanges(start = genes$start, end = genes$end),
  strand = ifelse(genes$strand == 1, "+", ifelse(genes$strand == -1, "-", "*")),
  feature_id = rownames(genes)
)
genes$chr <- genes$start <- genes$end <- genes$strand <- NULL

S4Vectors::elementMetadata(rowRanges)[names(genes)] <- genes

sce <- SummarizedExperiment(
  assays = list(counts = ge[["RNA"]]$counts),
  rowRanges = rowRanges,
  colData = ge@meta.data
)

pbge <- scsajr::pseudobulk(sce, c("sample_id", "celltype")) # calculate pseudobulk by celltype
SummarizedExperiment::assay(pbge, "cpm") <- scsajr::calc_cpm(pbge)
```

## 4.1 Visualise markers

``` r
markers <- scsajr::select_all_markers(pbasf@metadata$markers, pbasf@metadata$all_celltype_test, dpsi_thr = 0.5, n = 4)
# par(mfrow = c(1, 2), mar = c(3, 6, 1, 6), oma = c(0, 12, 0, 0))
par(mfrow = c(1, 2), mar = c(2, 4, 1, 4), oma = c(0, 4, 0, 0))
scsajr::marker_heatmap(
  pbasf, pbge, "celltype",
  psi_scale = FALSE,
  cpm_scale = TRUE,
  markers = markers,
  gene_names_col = "name"
)
```

    ## Warning in stats::cor(t(psi), use = "pairwise.complete.obs"): the standard
    ## deviation is zero

    ## Loading required package: randomcoloR

![](10x_gbm_files/figure-gfm/Visualise%20markers-1.png)<!-- -->

## 4.2 Coverage plots

\_from pre-saved files we can replot coverage for one of examples
generated by the pipeline

``` r
sid <- "ENSG00000092841.s17"
covs <- readRDS(paste0("./nf-scsajr_output/rds/examples_coverage/", sid, ".rds"))
scsajr::plot_segment_coverage(
  sid,
  chr = covs[[1]]$chr,
  start = covs[[1]]$start,
  stop = covs[[1]]$end,
  covs = covs,
  data_as = pbasf,
  data_ge = pbge,
  groupby = "celltype",
  barcodes = barcodes,
  samples = samples,
  gene_descr = gene_descr,
  plot_junc_only_within = FALSE,
  ylim_by_junc = TRUE,
  gtf = gtf,
  oma = c(4, 4, 4, 1)
)
```

![](10x_gbm_files/figure-gfm/Coverage%20plots-1.png)<!-- -->

We can zoom in on the alternative exon and focus on the most divergent
cell types.

``` r
scsajr::plot_segment_coverage(
  sid,
  chr = covs[[1]]$chr,
  start = 56160210,
  stop = 56161500,
  covs = covs,
  data_as = pbasf,
  data_ge = pbge,
  celltypes = c("6", "9", "5"),
  groupby = "celltype",
  barcodes = barcodes,
  samples = samples,
  gene_descr = gene_descr,
  plot_junc_only_within = FALSE,
  ylim_by_junc = TRUE,
  gtf = gtf,
  oma = c(4, 4, 4, 1)
)
```

![](10x_gbm_files/figure-gfm/Zoom%20in%20on%20alternative%20exon-1.png)<!-- -->

## 4.3 Differential splicing

Let’s compare just two celltypes, we will take cluster 5 and 6 as they
demonstrated divergent splicing in MYL6

``` r
cl5to6 <- scsajr::test_pair_as(pbasf, "celltype", c("5", "6"))
cl5to6[!is.na(cl5to6$dpsi) & cl5to6$fdr < 0.05 & abs(cl5to6$dpsi) > 0.9, ]
```

    ##                     overdispersion           pv          fdr       dpsi
    ## ENSG00000115935.s17             NA 2.176204e-15 6.772348e-13 -0.9705882
    ## ENSG00000080822.s18             NA 4.869223e-21 2.841191e-18  0.9312624
    ## ENSG00000217930.s9              NA 5.794861e-12 1.176105e-09 -1.0000000
    ## ENSG00000109083.s12             NA 1.849033e-11 3.319725e-09  0.9285714
    ## ENSG00000159267.s13             NA 8.734963e-09 8.321389e-07 -1.0000000

``` r
sid <- "ENSG00000217930.s9"
pbasf@metadata$all_celltype_test[sid, ]
```

    ##                    overdispersion        group   group_fdr low_state high_state
    ## ENSG00000217930.s9             NA 9.181843e-15 1.67692e-13         4          1
    ##                         dpsi
    ## ENSG00000217930.s9 0.3285871

``` r
pbasf@metadata$markers$fdr[sid, ]
```

    ##         1        10         2         3         4         5         6         7 
    ## 0.6862076 0.9610627 0.8694858 0.7374354 0.7678948 0.7883998 0.6406943 0.9618651 
    ##         8         9 
    ## 0.9289754 0.9712656

``` r
# _coverage plots from BAM file ####
# when segment of interest was not plotted by pipeline
# or if celltype annotation has changed, one may need to plot coverage from BAM files
covs <- NULL # plot_segment_coverage returns coverage object that can be used next time
covs <- scsajr::plot_segment_coverage(
  sid,
  start = 4343250,
  stop = 4351500,
  covs = covs,
  data_as = pbasf,
  data_ge = pbge,
  groupby = "celltype",
  barcodes = barcodes,
  samples = samples,
  gene_descr = gene_descr,
  celltypes = c("5", "6"),
  plot_junc_only_within = FALSE,
  ylim_by_junc = TRUE,
  gtf = gtf,
  oma = c(4, 4, 4, 1)
)
```

    ## Loading required package: GenomicAlignments

    ## Loading required package: Biostrings

    ## Loading required package: XVector

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

    ## Loading required package: Rsamtools

![](10x_gbm_files/figure-gfm/Inspect%20segment%20and%20plot-1.png)<!-- -->

``` r
# Exon seems to be included in cluster 5 and excluded in cluster 6. However, coverage is quite low
```

``` r
# We can save coverage for future use.
# It is not very important in case of single sample, but with tens of samples and celltypes extracting coverage from BAM files can take a while.
saveRDS(covs, paste0("./nf-scsajr_output/rds/examples_coverage/", sid, ".rds"))
```

## 4.4 domain enrichment analyses

``` r
# Load segment to domain annotation and domain description
domain2seg <- readRDS(url("https://github.com/cellgeni/nf-scsajr/raw/refs/heads/main/ref/human_2020A/functional_annotation/domain2seg.df.rds"))
domain_descr <- readRDS(url("https://github.com/cellgeni/nf-scsajr/raw/refs/heads/main/ref/human_2020A/functional_annotation/all_domain_descr.rds"))

# Load all segments
pbas_all <- readRDS("./nf-scsajr_output/rds/pbas.rds")

# Take all domains that have description
domain2seg <- domain2seg$all
domain2seg <- domain2seg[domain2seg$domain %in% domain_descr$ENTRY_AC, ]

domain_sites2use <- "ad" # focus on cassette exons only
# Select background set
# All genes that have at least one segment with reasonable coverage
gene_uni <- scsajr::filter_segments_and_samples(pbas_all, seg_min_sd = 0, celltype_min_samples = 1, sample_min_ncells = 30)
gene_uni <- unique(rowData(gene_uni)$gene_id)
seg_uni <- rownames(pbas_all)[segs$is_exon & segs$sites %in% domain_sites2use & segs$gene_id %in% gene_uni]

cl5to6sgn <- cl5to6[!is.na(cl5to6$dpsi) & cl5to6$fdr < 0.05 & abs(cl5to6$dpsi) > 0.1, ]
sids <- list(
  cl6 = rownames(cl5to6sgn)[cl5to6sgn$dpsi > 0.1],
  cl5 = rownames(cl5to6sgn)[cl5to6sgn$dpsi < -0.1]
)

ipro <- clusterProfiler::compareCluster(sids,
  fun = "enricher",
  universe = seg_uni,
  pAdjustMethod = "BH",
  TERM2GENE = domain2seg,
  TERM2NAME = domain_descr[, c("ENTRY_AC", "ENTRY_NAME")]
)

clusterProfiler::dotplot(ipro)
```

![](10x_gbm_files/figure-gfm/Domain%20enrichment%20analyses-1.png)<!-- -->

## 4.5 re-pseudobulk

``` r
# Using new cell clusters
# redo pbasf, pbge, barcodes
raw <- readRDS("./nf-scsajr_output/rds/sample.rds")
barcodes_new <- ge@meta.data[c("sample_id", "barcode", "seurat_clusters")]

raw$i <- raw$i[, rownames(barcodes_new)]
raw$e <- raw$e[, rownames(barcodes_new)]
raw$seg <- segs
pbas_new <- scsajr::make_summarized_experiment(raw, barcodes_new)
pbas_new <- scsajr::pseudobulk(pbas_new, "seurat_clusters")
pbasf_new <- scsajr::filter_segments_and_samples(
  pbas_new,
  groupby = "seurat_clusters",
  seg_min_sd = 0.1,
  celltype_min_samples = 1,
  sample_min_ncells = 30
)
# str(pbasf_new)
```
