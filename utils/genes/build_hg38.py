#!/usr/bin/env python3
"""
Builds the hg38 gene annotation files of a ProteinPaint server from public sources, keeping
the version of every transcript accession (NM_004456.5, ENST00000320356.7):

  <out>/refGene.hg38.gz(.tbi)        RefSeq track, from UCSC ncbiRefSeqCurated
  <out>/gencode.v<N>.hg38.gz(.tbi)   GENCODE track, from UCSC wgEncodeGencodeCompV<N> and PseudoGeneV<N>
  <out>/genes.hg38.db                gene db
  <out>/diff.md, diff.json           what changed against the previous gene db

Usage:
  pip install pysam
  python3 utils/genes/build_hg38.py --out <dir> [--gencode 50] [--old-db genes.hg38.db] [--cache <dir>]

Without --gencode, the latest release at UCSC is used. Without --old-db, the gene db that is
currently published at proteinpaint.stjude.org/ppSupport is downloaded. The previous db is
needed for the diff report and for the chrM rows, which no source table below has.

Rules, all checked to reproduce the 2025-10-19 gene db from the same sources:
- RefSeq rows: every row of ncbiRefSeqCurated; default isoform: member of ncbiRefSeqSelect
  (RefSeq Select / MANE); description: ncbiRefSeqLink.product
- GENCODE rows: every row of the comprehensive and pseudogene sets; default isoform and gene2canonicalisoform:
  the Ensembl_canonical tag; description: kgXref
- genealias: HGNC alias and previous symbols, and the ENSG accession of a gene
- codingGenes: genes with a coding RefSeq transcript; gene2coord: span of the default RefSeq
  isoform, or the widest one of a gene that has no default

Network: hgdownload.soe.ucsc.edu, storage.googleapis.com (HGNC), proteinpaint.stjude.org.
"""

import argparse, csv, datetime, gzip, json, os, re, sqlite3, sys, urllib.request
from email.utils import parsedate_to_datetime

UCSC = 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/'
HGNC = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
OLD_DB = 'https://proteinpaint.stjude.org/ppSupport/genes.hg38.db'
# genes that every build must have a default RefSeq and GENCODE isoform for
LANDMARKS = ['TP53', 'KRAS', 'EZH2', 'MYC', 'BRAF', 'EGFR', 'RUNX1', 'SETD2', 'KMT2D', 'BCR', 'ABL1']

source_dates = {}  # file name -> Last-Modified date of the download


def fetch(url, cache):
	"""download to the cache dir unless already there; returns the local path"""
	path = os.path.join(cache, url.rsplit('/', 1)[1])
	datefile = path + '.date'
	if not os.path.exists(path) or not os.path.exists(datefile):
		print('downloading', url, file=sys.stderr)
		with urllib.request.urlopen(url) as r, open(path + '.part', 'wb') as f:
			while chunk := r.read(1 << 20):
				f.write(chunk)
			lm = r.headers.get('Last-Modified')
		os.replace(path + '.part', path)
		with open(datefile, 'w') as f:
			f.write(parsedate_to_datetime(lm).date().isoformat() if lm else 'unknown')
	source_dates[os.path.basename(path)] = open(datefile).read().strip()
	return path


def rows(path):
	with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as f:
		for line in f:
			yield line.rstrip('\n').split('\t')


def latest_gencode():
	with urllib.request.urlopen(UCSC) as r:
		versions = re.findall(r'wgEncodeGencodeCompV(\d+)\.txt\.gz', r.read().decode())
	if not versions:
		raise SystemExit('no wgEncodeGencodeCompV* table found at UCSC')
	return max(map(int, versions))


def base(accession):
	return re.sub(r'\.\d+$', '', accession)


def gene_model(l, name, description=None, category=None):
	"""a row of a UCSC genePredExt table -> (chr, start, stop, model); same model as
	parseBedLine() of server/src/bedj.parseBed.js builds from a bigGenePred line"""
	isoform, chrom, strand = l[1], l[2], l[3]
	start, stop, cds_start, cds_stop = int(l[4]), int(l[5]), int(l[6]), int(l[7])
	starts = [int(i) for i in l[9].rstrip(',').split(',')]
	stops = [int(i) for i in l[10].rstrip(',').split(',')]
	exon, intron, coding, left, right = [], [], [], [], []
	for i, (a, b) in enumerate(zip(starts, stops)):
		exon.append([a, b])
		if i:
			intron.append([stops[i - 1], a])
		if a < cds_start:
			left.append([a, min(b, cds_start)])
		if max(a, cds_start) < min(b, cds_stop):
			coding.append([max(a, cds_start), min(b, cds_stop)])
		if b > cds_stop:
			right.append([max(a, cds_stop), b])
	forward = strand == '+'
	if not forward:  # every list goes 5' to 3'
		for lst in (exon, intron, coding, left, right):
			lst.reverse()
	m = {'name': name, 'isoform': isoform, 'strand': strand, 'exon': exon, 'rnalen': sum(b - a for a, b in exon)}
	if category:
		m['category'] = category
	if description:
		m['description'] = description
	if intron:
		m['intron'] = intron
	if cds_start != cds_stop:
		utr5, utr3 = (left, right) if forward else (right, left)
		m.update(cdslen=sum(b - a for a, b in coding), codingstart=cds_start, codingstop=cds_stop, coding=coding)
		if utr5:
			m['utr5'] = utr5
		if utr3:
			m['utr3'] = utr3
		# as checkReadingFrame.js: a first coding exon that does not begin a codon
		frames = l[15].rstrip(',').split(',') if len(l) > 15 else []
		pos = cds_start if forward else cds_stop
		for i, (a, b) in enumerate(zip(starts, stops) if forward else reversed(list(zip(starts, stops)))):
			if a <= pos <= b:
				idx = i if forward else len(starts) - 1 - i
				if idx < len(frames) and frames[idx] in ('1', '2'):
					m['startCodonFrame'] = int(frames[idx])
				break
	return chrom, start, stop, m


def write_track(items, path):
	import pysam

	items.sort(key=lambda i: (i[0], i[1], i[2]))
	plain = path[:-3]
	with open(plain, 'w') as f:
		for chrom, start, stop, m in items:
			f.write(f'{chrom}\t{start}\t{stop}\t{json.dumps(m, separators=(",", ":"))}\n')
	pysam.tabix_compress(plain, path, force=True)
	os.remove(plain)
	pysam.tabix_index(path, preset='bed', force=True)


def is_primary(chrom):
	return '_' not in chrom


def build(args):
	os.makedirs(args.out, exist_ok=True)
	os.makedirs(args.cache, exist_ok=True)
	v = args.gencode or latest_gencode()
	f = lambda name: fetch(UCSC + name + '.txt.gz', args.cache)
	old_db = args.old_db or fetch(OLD_DB, args.cache)

	# ---------- RefSeq
	product = {l[0]: l[3] for l in rows(f('ncbiRefSeqLink'))}
	select = {l[1] for l in rows(f('ncbiRefSeqSelect'))}
	refseq = [gene_model(l, l[12], product.get(l[1])) for l in rows(f('ncbiRefSeqCurated'))]

	# ---------- GENCODE
	attrs = {l[4]: l for l in rows(f(f'wgEncodeGencodeAttrsV{v}'))}  # transcriptId -> row
	canonical = {l[0] for l in rows(f(f'wgEncodeGencodeTagV{v}')) if l[1] == 'Ensembl_canonical'}
	kgdesc = {l[0]: l[7] for l in rows(f('kgXref'))}
	gencode = []
	for l in [*rows(f(f'wgEncodeGencodeCompV{v}')), *rows(f(f'wgEncodeGencodePseudoGeneV{v}'))]:
		a = attrs.get(l[1])
		desc, category = kgdesc.get(l[1]), a[12] if a else None
		is_canonical = l[1] in canonical
		# the chrY copy of a PAR transcript is the same accession, as in RefSeq
		l[1] = l[1].replace('_PAR_Y', '')
		gencode.append(gene_model(l, l[12], desc, category) + (is_canonical,))

	write_track(list(refseq), os.path.join(args.out, 'refGene.hg38.gz'))
	write_track([g[:4] for g in gencode], os.path.join(args.out, f'gencode.v{v}.hg38.gz'))

	# ---------- gene db
	dbfile = os.path.join(args.out, 'genes.hg38.db')
	if os.path.exists(dbfile):
		os.remove(dbfile)
	db = sqlite3.connect(dbfile)
	db.executescript(
		"""
CREATE TABLE genes (name character varying(50) collate nocase, isoform character varying(50) collate nocase, isdefault smallint, genemodel json);
CREATE TABLE genealias (alias character varying(50) collate nocase, name character varying(50) collate nocase);
CREATE TABLE refseq2ensembl (ensembl char(100) collate nocase, refseq char(100) collate nocase);
CREATE TABLE gene2coord (name char(100) collate nocase, chr char(100) collate nocase, start INTEGER, stop INTEGER);
CREATE TABLE gene2canonicalisoform (gene character collate nocase, isoform character collate nocase);
CREATE TABLE codingGenes (name character collate nocase);
CREATE TABLE buildDate (date char(100), sources text);
"""
	)

	def db_row(chrom, start, stop, m, isdefault, trackname):
		return m['name'], m['isoform'], int(isdefault), json.dumps({**m, 'chr': chrom, 'start': start, 'stop': stop, 'trackname': trackname}, separators=(',', ':'))

	genes = [db_row(*g, g[3]['isoform'] in select, 'RefGene') for g in refseq]
	genes += [db_row(*g[:4], g[4], 'GENCODE') for g in gencode]
	# the one row of the previous db that no source table has: the chrM polycistronic transcript
	old = sqlite3.connect(f'file:{old_db}?mode=ro', uri=True)
	genes += old.execute("select name,isoform,isdefault,genemodel from genes where isoform='mt'").fetchall()
	db.executemany('insert into genes values (?,?,?,?)', genes)

	alias = set()
	with open(fetch(HGNC, args.cache), encoding='utf-8') as fh:
		for r in csv.DictReader(fh, delimiter='\t'):
			for key in ('alias_symbol', 'prev_symbol'):
				alias.update((a, r['symbol']) for a in r[key].strip('"').split('|') if a)
			if r['ensembl_gene_id']:
				alias.add((r['ensembl_gene_id'], r['symbol']))
	alias.update((base(a[0]).replace('_PAR_Y', ''), a[1]) for a in attrs.values())  # ENSG -> symbol
	db.executemany('insert into genealias values (?,?)', sorted(alias))

	# ENST unversioned, as the consumer queries the gdc api with it; refseq in both forms
	r2e = set()
	for l in rows(f(f'wgEncodeGencodeRefSeqV{v}')):
		r2e.update({(base(l[0]), l[1]), (base(l[0]), base(l[1]))})
	db.executemany('insert into refseq2ensembl values (?,?)', sorted(r2e))

	coord = {}  # gene -> (rank, chr, start, stop); the default isoform on a primary chr wins
	for chrom, start, stop, m in refseq:
		rank = (m['isoform'] not in select, not is_primary(chrom), start - stop)
		if m['name'] not in coord or rank < coord[m['name']][0]:
			coord[m['name']] = (rank, chrom, start, stop)
	db.executemany('insert into gene2coord values (?,?,?,?)', [(k, *c[1:]) for k, c in coord.items()])

	g2c = {}
	for chrom, start, stop, m, is_canonical in gencode:
		if is_canonical and is_primary(chrom):
			g2c.setdefault(m['name'], m['isoform'])
			a = attrs.get(m['isoform'])
			if a:
				g2c.setdefault(base(a[0]), m['isoform'])
	db.executemany('insert into gene2canonicalisoform values (?,?)', g2c.items())

	coding = sorted({m['name'] for _, _, _, m in refseq if 'cdslen' in m})
	db.executemany('insert into codingGenes values (?)', [(n,) for n in coding])

	today = datetime.date.today().isoformat()
	sources = f"RefSeq: UCSC ncbiRefSeqCurated {source_dates['ncbiRefSeqCurated.txt.gz']}; GENCODE v{v}; HGNC {source_dates['hgnc_complete_set.txt']}"
	db.execute('insert into buildDate values (?,?)', (today, sources))
	db.executescript(
		"""
CREATE INDEX genesymbol on genes(name collate nocase);
CREATE INDEX geneisoform on genes(isoform collate nocase);
CREATE INDEX aliastext on genealias(alias collate nocase);
CREATE INDEX r2er on refseq2ensembl(refseq);
CREATE INDEX r2ee on refseq2ensembl(ensembl);
CREATE INDEX gname on gene2coord(name);
CREATE INDEX gene2canonicalisoform_gene on gene2canonicalisoform(gene);
"""
	)
	db.commit()

	selfcheck(db)
	diff(old, db, args.out, today, sources)
	print(f'built {dbfile}: {sources}', file=sys.stderr)


def selfcheck(db):
	"""stop the build rather than publish a db that would break a server"""
	q = lambda sql, *p: db.execute(sql, p).fetchone()[0]
	assert q("select count(*) from genes where isoform like 'N%'") > 90000, 'too few RefSeq rows'
	assert q("select count(*) from genes where isoform like 'ENST%'") > 300000, 'too few GENCODE rows'
	assert q("select count(*) from genes where isoform like 'N%' and isoform not like '%.%'") == 0, 'RefSeq accession without version'
	assert q("select count(*) from genes where isoform like 'ENST%' and isoform not like '%.%'") == 0, 'ENST accession without version'
	assert q('select count(*) from genealias') > 100000, 'too few aliases'
	assert q('select count(*) from codingGenes') > 19000, 'too few coding genes'
	for gene in LANDMARKS:
		for prefix in ('NM_%', 'ENST%'):
			n = q('select count(*) from genes where name=? and isdefault=1 and isoform like ?', gene, prefix)
			assert n >= 1, f'{gene} has no default {prefix[:-1]} isoform'
		assert q('select count(*) from gene2coord where name=?', gene) == 1, f'{gene} missing from gene2coord'
		assert q('select count(*) from gene2canonicalisoform where gene=?', gene) == 1, f'{gene} has no canonical isoform'


def diff(old, new, out, today, sources):
	"""compares the two dbs transcript by transcript, matched on the accession without version"""

	def load(db):
		d = {}
		for name, isoform, isdefault, gm in db.execute('select name,isoform,isdefault,genemodel from genes'):
			m = json.loads(gm)
			if is_primary(m.get('chr', '')):  # an alt contig copy says nothing new
				d.setdefault(base(isoform), (name, isoform, isdefault, m))
		return d

	a, b = load(old), load(new)
	shared = [k for k in a if k in b]
	changed = lambda key: [k for k in shared if key(a[k]) != key(b[k])]
	model = changed(lambda r: (r[3]['exon'], r[3].get('coding')))
	protein = changed(lambda r: r[3].get('coding'))
	aa = lambda m: f'{m["cdslen"] // 3 - 1} aa' if 'cdslen' in m else 'noncoding'
	lists = {
		'added': sorted(b[k][1] for k in b if k not in a),
		'removed': sorted(a[k][1] for k in a if k not in b),
		# empty the first time a versioned db replaces an unversioned one
		'new_version': sorted(f'{a[k][1]} -> {b[k][1]}' for k in shared if '.' in a[k][1] and a[k][1] != b[k][1]),
		'model_changed': sorted(b[k][1] for k in model),
		'protein_changed': sorted(f'{b[k][0]} {b[k][1]}: {aa(a[k][3])} -> {aa(b[k][3])}' for k in protein),
		'symbol_renamed': sorted(f'{a[k][0]} -> {b[k][0]} ({b[k][1]})' for k in changed(lambda r: r[0].upper())),
		'default_changed': sorted(f'{b[k][0]} {b[k][1]}: {"gained" if b[k][2] else "lost"} default' for k in changed(lambda r: bool(r[2]))),
	}
	old_date = old.execute('select date from buildDate').fetchone()[0]
	summary = {'built': today, 'sources': sources, 'previous_build': old_date, 'counts': {k: len(v) for k, v in lists.items()}}
	with open(os.path.join(out, 'diff.json'), 'w') as fh:
		json.dump({**summary, **lists}, fh, indent=1)
	labels = {
		'added': 'Transcripts added',
		'removed': 'Transcripts removed',
		'new_version': 'Transcripts with a new version',
		'model_changed': 'Transcripts whose exons or CDS moved',
		'protein_changed': '... of which the CDS changed',
		'symbol_renamed': 'Transcripts whose gene symbol changed',
		'default_changed': 'Transcripts that gained or lost default status',
	}
	with open(os.path.join(out, 'diff.md'), 'w') as fh:
		fh.write(f'# hg38 gene models {today}\n\n{sources}\n\nCompared with the gene db built {old_date}, primary chromosomes only.\n\n| Change | Count |\n|---|---:|\n')
		fh.writelines(f'| {labels[k]} | {len(v)} |\n' for k, v in lists.items())
		for k in ('protein_changed', 'symbol_renamed', 'default_changed'):
			fh.write(f'\n## {labels[k].lstrip(". ")} (first 50; all in diff.json)\n\n')
			fh.writelines(f'- {i}\n' for i in lists[k][:50])


if __name__ == '__main__':
	p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	p.add_argument('--out', required=True, help='output dir')
	p.add_argument('--gencode', type=int, help='GENCODE release, default: latest at UCSC')
	p.add_argument('--old-db', help='previous genes.hg38.db, default: download the published one')
	p.add_argument('--cache', help='dir for downloads, default: <out>/downloads')
	args = p.parse_args()
	args.cache = args.cache or os.path.join(args.out, 'downloads')
	build(args)
