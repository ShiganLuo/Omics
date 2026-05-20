1. StringTie's primary GTF output
The primary output of StringTie is a Gene Transfer Format (GTF) file that contains details of the transcripts that StringTie assembles from RNA-Seq data. GTF is an extension of GFF (Gene Finding Format, also called General Feature Format), and is very similar to GFF2 and GFF3. The field definitions for the 9 columns of GTF output can be found at the Ensembl site here. The following is an example of a transcript assembled by StringTie as shown in a GTF file (scroll right within the box to see the full field contents):
Description of each column's values:

```plain
seqname source      feature     start   end     score   strand  frame attributes
chrX    StringTie   transcript  281394  303355  1000    +       .     gene_id "ERR188044.1"; transcript_id "ERR188044.1.1"; reference_id "NM_018390"; ref_gene_id "NM_018390"; ref_gene_name "PLCXD1"; cov "101.256691"; FPKM "530.078918"; TPM "705.667908";
chrX    StringTie   exon        281394  281684  1000    +       .     gene_id "ERR188044.1"; transcript_id "ERR188044.1.1"; exon_number "1"; reference_id "NM_018390"; ref_gene_id "NM_018390"; ref_gene_name "PLCXD1"; cov "116.270836";
...
```
- seqname: Denotes the chromosome, contig, or scaffold for this transcript. Here the assembled transcript is on chromosome X.
- source: The source of the GTF file. Since this example was produced by StringTie, this column simply shows 'StringTie'.
- feature: Feature type; e.g., exon, transcript, mRNA, 5'UTR.
- start: Start position of the feature (exon, transcript, etc), using a 1-based index.
- end: End position of the feature, using a 1-based index.
- score: A confidence score for the assembled transcript. Currently this field is not used, and StringTie reports a constant value of 1000 if the transcript has a connection to a read alignment bundle.
- strand: If the transcript resides on the forward strand, '+'. If the transcript resides on the reverse strand, '-'.
- frame: Frame or phase of CDS features. StringTie does not use this field and simply records a ".".
- attributes: A semicolon-separated list of tag-value pairs, providing additional information about each feature. Depending on whether an instance is a transcript or an exon and on whether the transcript matches the reference annotation file provided by the user, the content of the attributes field will differ. The following list describes the possible attributes shown in this column:
    - gene_id: A unique identifier for a single gene and its child transcript and exons based on the alignments' file name.
    - transcript_id: A unique identifier for a single transcript and its child exons based on the alignments' file name.
    - exon_number: A unique identifier for a single exon, starting from 1, within a given transcript.
    - reference_id: The transcript_id in the reference annotation (optional) that the instance matched.
    - ref_gene_id: The gene_id in the reference annotation (optional) that the instance matched.
    - ref_gene_name: The gene_name in the reference annotation (optional) that the instance matched.
    - cov: The average per-base coverage for the transcript or exon.
    - FPKM: Fragments per kilobase of transcript per million read pairs. This is the number of pairs of reads aligning to this feature, normalized by the total number of fragments sequenced (in millions) and the length of the transcript (in kilobases).
    - TPM: Transcripts per million. This is the number of transcripts from this particular gene normalized first by gene length, and then by sequencing depth (in millions) in the sample. A detailed explanation and a comparison of TPM and FPKM can be found here, and TPM was defined by B. Li and C. Dewey here.

2. Gene abundances in tab-delimited format
If StringTie is run with the -A <gene_abund.tab> option, it returns a file containing gene abundances. The tab-delimited gene abundances output file has nine fields; here is an example of a gene abundance file produced by StringTie using reference annotation:

```plain
Gene ID     Gene Name   Reference   Strand  Start   End     Coverage    FPKM        TPM
NM_000451   SHOX        chrX        +       624344  646823  0.000000    0.000000    0.000000
NM_006883   SHOX        chrX        +       624344  659411  0.000000    0.000000    0.000000
...
```
- Column 1 / Gene ID: The gene identifier comes from the reference annotation provided with the -G option. If no reference is provided this field is replaced with the name prefix for output transcripts (-l).
- Column 2 / Gene Name: This field contains the gene name in the reference annotation provided with the -G option. If no reference is provided this field is populated with '-'.
- Column 3 / Reference: Name of the reference sequence that was used in the alignment of the reads. Equivalent to the 3rd column in the .SAM alignment.
- Column 4 / Strand: '+' denotes that the gene is on the forward strand, '-' for the reverse strand.
- Column 5 / Start: Start position of the gene (1-based index).
- Column 6 / End: End position of the gene (1-based index).
- Column 7 / Coverage: Per-base coverage of the gene.
- Column 8 / FPKM: normalized expression level in FPKM units (see previous section).
- Column 9 / TPM: normalized expression level in RPM units (see previous section).

3. Fully covered transcripts matching the reference annotation transcripts (in GTF format)

If StringTie is run with the -C <cov_refs.gtf> option (requires -G <reference_annotation>), it returns a file with all the transcripts in the reference annotation that are fully covered, end to end, by reads. The output format is a GTF file as described above. Each line of the GTF is corresponds to a gene or transcript in the reference annotation.

4. Ballgown Input Table Files

If StringTie is run with the -B option, it returns a Ballgown input table file, which contains coverage data for all transcripts. The output table files are placed in the same directory as the main GTF output. These tables have these specific names: (1) e2t.ctab, (2) e_data.ctab, (3) i2t.ctab, (4) i_data.ctab, and (5) t_data.ctab. A detailed description of each of these five required inputs to Ballgown can be found on the Ballgown site, at this link.

5. Merge mode: Merged GTF

If StringTie is run with the --merge option, it takes as input a list of GTF/GFF files and merges/assembles these transcripts into a non-redundant set of transcripts. This step creates **a uniform set of transcripts for all samples** to facilitate the downstream calculation of differentially expressed levels for all transcripts among the different experimental conditions. Output is a merged GTF file with all merged gene models, but without any numeric results on coverage, FPKM, and TPM. Then, with this merged GTF, StringTie can re-estimate abundances by running it again with the -e option on the original set of alignment files, as illustrated in the figure below.
