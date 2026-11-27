---
name: update-gene-models
description: Rebuild the hg38 gene annotation of ProteinPaint (RefSeq, GENCODE, HGNC aliases) from the latest public releases, with transcript versions kept, and review what changed. Use for the yearly gene model update, when asked to update GENCODE or RefSeq, or when the "Update gene models" workflow has opened an issue.
---

# Update the hg38 gene models

Done once a year, in the fall: NCBI publishes one GRCh38 RefSeq annotation release a year (August),
UCSC picks it up a few weeks later, and GENCODE has released twice by then.

The build is `utils/genes/build_genes.py`. Read its docstring first: it lists the sources and the
rule behind each table. It writes `refGene.hg38.gz`, `gencode.v<N>.hg38.gz`, `genes.hg38.db`,
`proteindomain.db`, `diff.md` and `diff.json`. Accessions keep their version (`NM_004456.5`); the server finds a
transcript at any version or none (`resolveIsoform()` in `server/src/genedbLookups.ts`).

## Steps

1. Build, either here or by running the "Update gene models" workflow and downloading its artifact:
   ```
   pip install pysam
   python3 utils/genes/build_genes.py --genome hg38 --out <dir outside the repo>
   ```
   It needs `hgdownload.soe.ucsc.edu`, `ftp.ncbi.nlm.nih.gov`, `storage.googleapis.com` and
   `proteinpaint.stjude.org`, about
   3 GB of disk and 6 minutes. It stops on its own if the result would break a server (`selfcheck()`).
   If a host is blocked, say which one and stop; do not build from a partial set of sources.

2. Read `diff.md` and judge it, do not just pass it on. Expect a few hundred to a few thousand
   moved models from a GENCODE bump and few or none from RefSeq. Look into, and report:
   - any gene in "the CDS changed" that is a known cancer gene: its residue numbers move for
     every dataset annotated on the previous release
   - "gene symbol changed": renamed genes drop out of gene sets and expression matrices
     that still use the old symbol, unless the old one is now in `genealias`
   - "Transcripts removed" far above the previous year's count, which means a source table moved
     or changed format rather than that biology changed
   - a build stopped by "the MANE release does not match this annotation": MANE names exact
     versions, so its current release must be the one made from the RefSeq and GENCODE releases
     being built. Report which of the three is behind; do not loosen the check
   - a RefSeq source date that has not moved since the previous build: then only GENCODE and
     HGNC are new, and the report should say so

3. Run the `check-transcript-versions` skill on the previous and the new gene db.

4. Hand over. You cannot publish the files: a person uploads the gene db, the two tracks and
   their `.tbi` to `ppSupport` and to each server's `tp/anno/`, and `proteindomain.db` to
   `tp/anno/db/`. The gene db and the domain db go together: the domain db is keyed by the
   versioned accessions of the gene db it was built with. Give them `diff.md` and your findings.

5. Only after they confirm the upload, open a PR that:
   - sets the GENCODE file name in `server/genome/hg38.base.ts` and
     `utils/public-support-files-download.sh` to the new release (a server whose `tp/anno/`
     lacks the file does not start, which is why this waits for the upload)
   - adds a line to `release.txt` naming the RefSeq source date and the GENCODE release

## hg19 and mm10

The same script builds them (`--genome hg19`, `--genome mm10`), but their sources are frozen:
rebuild them only when the script changes, not every year. hg19 shares the protein domain db with
hg38, so build hg38 first and pass `--domains <hg38 out>/proteindomain.db` to the hg19 build.
mm10 has no MANE, no HGNC aliases, no protein domains, and no default GENCODE isoform, since
GENCODE M25 does not mark a canonical transcript. After the upload, hg19 needs its GENCODE file
name in `server/genome/hg19.ts`, and mm10 a GENCODE track added to `server/genome/mm10.ts`.

## Not covered

- A GENCODE transcript with no RefSeq transcript of the same coding exons has no protein domains,
  nor does a retired version, nor a transcript whose protein NCBI has revised since the RefSeq
  release being built. The report line "Protein domains" gives the counts; a drop of more than a
  few hundred RefSeq transcripts against the previous year means NCBI's weekly protein files have
  moved well ahead of the release UCSC serves, and is worth reporting.

- The GDC gene db (GENCODE v36, frozen with each GDC data release) and `genes.hg38.mmrf.db`
  (GENCODE v32) are pinned on purpose. Do not rebuild them here.
- mm9, dm3, dm6, rn6 and the other genomes have no builder in this repo.
- Retired RefSeq versions (`ncbiRefSeqHistorical`) go in the `retiredgenes` table and open only
  when asked for by their exact accession. GENCODE has no such table: a retired ENST version
  opens the current model, with a note in the gene view.
- The test fixtures under `server/test/tp/anno/` stay as they are: they are the unversioned
  layout that the tests of the old form run against.
