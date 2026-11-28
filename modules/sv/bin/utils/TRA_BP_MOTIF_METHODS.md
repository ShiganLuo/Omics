# TRA Breakpoint Motif Mining — Methods (stage 4-7)

Companion to `tra_bp_motif.py` + `tra_bp_motif_figures.py`; continues the
`tra_te_analysis.py` evidence funnel (T0-T3) with a breakpoint-sequence layer:
background composition, motif enrichment, junction microhomology, and per-event
mechanism calls.

## Inputs

| input | path (as run) | role |
|---|---|---|
| events | `.../youth_vs_aged/tra_te/tra_events_annotated.tsv` | 143 TRA events + TE annotation (from `sv_tra_te.sh`) |
| homology | `.../youth_vs_aged/tra_te/stage3_homology.tsv` | NAHR tier (blastn ident/HSP per event) |
| reference | `/data/16T/Database/Reference/mulatta/NCBI/Mmul_T2T/GCF_049350105.2_T2T-MMU8v2.0_genomic.fa` (+`.fai`) | window extraction (fai random access, no pysam) |
| TE BED | `.../Mmul_T2T/rheMac10_rmsk_TE.lifted.bed` (cols 1-4,7,8) | TE-context matching + TE class labels |
| telocen BED | `.../Mmul_T2T.telocen_regions.bed` | telomere/centromere context |
| sample VCFs | `variation/germline_sv/<sample>/<sample>.sv.vcf.gz` | recover BND join form (merged VCF lost STRANDS) |
| motif lib | `bp_motif/motif_lib/*.motif` | HOMER-format PWMs (ctcf, ctcf-mys, yy1, sp1, nrf1, e2f1, ap1, cmyc, ets1, usf1), snapshotted from `/data/16T/Database/env/homer/homer.sif` |

## Design

1. **Breakpoint clusters.** 286 breakpoints of 143 events are clustered within
   500 bp per chromosome (246 clusters) so near-duplicate calls (SURVIVOR
   CIPOS ±109) do not pseudo-replicate the enrichment tests. Window = ±250 bp
   (501 bp) around each cluster representative; windows <60% of full length
   (chromosome ends) are dropped.
2. **Two backgrounds**, K=10 windows per representative, positions >=10 kb from
   any breakpoint:
   - `bg_uniform`: uniform on the same chromosome (total enrichment view);
   - `bg_matched`: same TE-context (in-TE vs non-TE at the exact position) and
     GC within ±5 pp of the breakpoint window (residual enrichment view,
     controls for GC and generic TE membership — NOT for TE family; see the
     strata check below).
3. **Composition** (stage4): GC, CpG O/E, max homopolymer, max mono/di/tri/tetra
   tandem run, purine-tract length, max AT fraction in any 30 bp —
   Mann-Whitney U vs each background.
4. **Structural motifs** (stage5a), presence per window:
   - G4: `(G{3,}N{1,7}){3}G{3,}` and the C-rich mirror;
   - L1 endonuclease consensus `TTTT/AA` — strict `TTTTAA|TTAAAA`, relaxed
     `TTTTA|TAAAA`;
   - Z-DNA: `(CG){4,}` or `(CA){5,}|(TG){5,}`;
   - inverted repeats: exact 8-mer seed whose downstream partner is its reverse
     complement with 0-12 bp spacer, exact outward extension (stem >=8 bp);
   - direct repeats: same seed with 1-20 bp spacer;
   - AG-mirror `(AG|GA){6,}` (H-DNA substrate), poly-A/T `A{10,}|T{10,}`,
     telomere repeats `TTAGGG|CCCTAA` (single and x2 tandem).
5. **TF PWMs** (stage5b): HOMER `.motif` matrices scanned both strands,
   log-odds vs uniform background, hit = best score >= 0.8 x max_score (the
   `motif_scan.py` convention).
6. **k-mers** (stage5c): canonical (strand-merged) 6-mer occurrence counts,
   Fisher exact on counts vs matched background, BH-FDR.
7. **Junction microhomology** (stage7): 40 bp reference flanks x/y (left/right
   of bp1) and u/v (bp2); the four BND join geometries
   A `t[p[` = x|v, B `t]p]` = x|rc(u), C `]p]t` = u|y, D `[p[t` = rc(v)|y;
   per-geometry longest suffix-prefix overlap (<=30 bp). The called geometry is
   taken from the per-sample BND ALT (matched within 200 bp on both mates);
   records without a sample match default to A and are flagged
   `merged_default`. Null: one matched-background position pair per event
   (same two clusters' matched sets), max over the four geometries.
8. **Mechanism call** (priority order): `TE-NAHR` if stage-3 blastn ident >=90%
   and HSP >=100 bp; else `microhomology-altEJ` if called-join MH >=4 bp; else
   `short-MH` (2-3 bp); else `blunt-NHEJ-like` (<=1 bp). Structural tags (G4 /
   Z-DNA / AG-mirror / telomere repeat / L1 EN at either breakpoint) are
   recorded separately; they are not used in the call.
9. **Stratified checks**: motif presence in NAHR vs non-NAHR clusters
   (`stage5_stratified.tsv`), and — via `windows_all.tsv` (per-window features +
   TE class + TF hits) — within TE-class strata (confound check for family
   composition).

## Statistics

Window-level presence, Fisher exact two-sided, odds ratio with Woolf 95% CI,
BH-FDR over all motifs per background; composition via Mann-Whitney U;
k-mers via count-level Fisher. All p-values are exploratory (n=143 events).

## Key numbers (run 2026-11-28, seed 42)

- 143 events -> 286 breakpoints -> 246 clusters; backgrounds 2304 matched /
  2460 uniform windows.
- Composition: CpG O/E 0.355 vs 0.236 matched (p=5.2e-21); GC 0.436 vs 0.426
  (p=0.039); homopolymer/tandem/purine/AT-rich: ns.
- Motifs (matched bg): G4 C-rich OR 3.70 (FDR 3.7e-4); G4 G-rich OR 3.26
  (FDR 2.6e-3); E2F1 PWM OR 9.65 (FDR 4.6e-4). L1 EN (strict and relaxed),
  Z-DNA, inverted/direct repeats, poly-A/T, telomere repeats, CTCF, YY1, Sp1,
  NRF1, c-Myc, AP-1, ETS1, USF1: not significant.
- TE-stratified confound check (`windows_all.tsv`): G4 G-rich survives in
  non-TE windows (OR 4.1, p=9.4e-4); G4 C-rich survives in LTR (OR 7.0,
  p=0.044) and non-TE (OR 3.8, p=9.7e-4); E2F1 survives in LTR (5/62 vs 0/278,
  p=1.8e-4) and non-TE (OR 7.1, p=0.024); CpG O/E elevated in every stratum.
  => signals are not ERV-LTR family-composition artifacts.
- Junction MH >=4 bp: 2/143 observed vs 0/143 null (Fisher p=0.50) — bulk
  repair is homology-free.
- Mechanism: TE-NAHR 19, microhomology-altEJ 2, short-MH 2, blunt-NHEJ-like
  120. Join forms (recovered from sample BNDs): t[p[ 71, t]p] 72.
- HOMER de novo + known enrichment (bg = matched windows, `homer/`): E-box
  (E2A/NGA4/Ptf1a/ZEB1), GC-rich zinc-finger (GLIS3/KLF1/ZNF416), TATA-box,
  poly-T — promoter/origin-like signature, consistent with the PWM layer.

## Rerun

One command: `bash .../youth_vs_aged/sv_tra_bp_motif.sh` (engine + HOMER de
novo + figures + `verify_te_strata.py`; HOMER step needs
`apptainer exec --bind /data/16T`). Review copies of all three scripts live in
`<outdir>/scripts/` (user rule: reproducible scripts ship with the output);
the repo files here are the canonical maintained versions. Individual stages
are the commands recorded in that script. Interpreter:
`/home/luosg/miniforge3/envs/plot/bin/python` (pandas/numpy/scipy/matplotlib).

## Limitations

- Junction sequence is NOT assembled (pbsv BND records carry no SVINSSEQ);
  microhomology is a reference-flank proxy for the junction overlap. Definitive
  junction evidence needs per-event read-level assembly (HiFi/fiber).
- `bg_matched` controls GC + TE membership but not TE family; the
  `windows_all.tsv` strata check addresses this, but family-level residual
  confounding cannot be fully excluded with n=246 windows.
- TF PWMs come from human/mouse HOMER matrices; a hit is a sequence-similarity
  call, not binding evidence. E2F1 hits are 8 windows — candidate level.
- 1v1 design (one young, one aged individual): group labels indicate presence
  per individual, not group biology; "shared" events may include germline
  polymorphism or recurrent caller artifacts.
- The called-join form is recovered for 143/143 events from sample BNDs within
  200 bp, but coordinates between caller and SURVIVOR merge can differ by up to
  the CIPOS window.
