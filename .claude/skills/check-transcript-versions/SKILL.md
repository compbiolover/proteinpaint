---
name: check-transcript-versions
description: Verify that ProteinPaint finds well known transcripts (TP53, EZH2, KRAS, MYC and others) whether they are asked for without a version, at the current version or at a retired one, against a gene db stored with or without versions. Use before a new gene db is deployed, after changing gene or isoform lookup code, or when a variant is reported at the wrong residue or an isoform is reported as not found.
---

# Check transcript version handling

A gene db stores isoform accessions either without a version (dbs built through 2025:
`NM_004456`) or with one (`utils/genes/build_hg38.py`: `NM_004456.5`). Requests, urls, saved
sessions and data files carry either form. All of them must reach the same transcript.

## Steps

1. Unit tests, which build both kinds of db in memory:
   ```
   cd server
   npx tsx --conditions=sjpp/dev src/test/genedbLookups.unit.spec.ts
   npx tsx --conditions=sjpp/dev src/test/geneRefValidation.unit.spec.ts
   cd ../shared/utils && npx tsx src/test/helpers.unit.spec.ts
   ```
   They need a `serverconfig.json`; `container/ci/serverconfig.json` with its paths pointed at this
   checkout works.

2. Real gene dbs, old against new. The published one is at
   `https://proteinpaint.stjude.org/ppSupport/genes.hg38.db`.
   ```
   cd server
   npx tsx --conditions=sjpp/dev ../utils/genes/checkLookups.ts <old genes.hg38.db> <new genes.hg38.db>
   ```
   Every line must be `ok`. A line marked "protein length differs" is not a failure of the
   lookup. When the accession on that line is a retired version, the new db opened the retired
   model, as it should (the EZH2 hotspot is Y641 on `NM_004456.1`, 746 aa, and Y646 on
   `NM_004456.5`, 751 aa). When it is a current or unversioned accession, the new release changed
   that protein: report it with both lengths, since every residue number downstream of the change
   moves for data annotated on the old release.

3. When adding a gene to the check, add it to `cases` in `utils/genes/checkLookups.ts` with one
   accession of each form.

## Where versions are handled

- The gene view shows a note (`mayShowVersionNote()` in `client/src/block.init.js`) when the
  transcript opened is a retired version, or is not the version that was asked for.
- A retired version has no protein domains: the domain db is positioned on the current version.

- `resolveIsoform()` in `server/src/genedbLookups.ts`: what was asked for -> what the db stores.
  Every isoform getter of `genome.genedb` goes through it.
- `sameAccession()` in `shared/utils/src/helpers.ts`: for comparing an accession of a data file
  or url with one of a gene model. Use it instead of `==` wherever the two sides can come from
  different releases.
- A dataset's own `byisoform` getter receives the accession as the gene db stores it. A getter
  that queries an external api or file keyed without versions must strip it with
  `stripAccessionVersion()`.
