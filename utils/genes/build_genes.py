#!/usr/bin/env python3
"""
Builds the gene annotation files of a ProteinPaint genome from public sources, keeping the
version of every transcript accession (NM_004456.5, ENST00000320356.7):

  <out>/refGene.<genome>.gz(.tbi)          RefSeq track, from UCSC ncbiRefSeqCurated
  <out>/gencode.v<N>.<genome>.gz(.tbi)     GENCODE track, from UCSC wgEncodeGencodeComp and PseudoGene
  <out>/genes.<genome>.db                  gene db
  <out>/proteindomain.db                   protein domains and sites of each coding transcript (human)
  <out>/diff.md, diff.json                 what changed against the previous gene db

Usage:
  pip install pysam
  python3 utils/genes/build_genes.py --genome hg38|hg19|mm10 --out <dir> [--gencode 50] [--old-db <file>] [--cache <dir>]

Without --gencode, the latest release that UCSC has for the genome is used (hg19: the latest lifted
to GRCh37). Without --old-db, the gene db that is currently published at
proteinpaint.stjude.org/ppSupport is downloaded: it is needed for the diff report and for the
few rows that no source table has. hg38 is the genome with new releases every year; the sources
of hg19 and mm10 are frozen, so those two are rebuilt only when this script changes.

hg38 and hg19 share one protein domain db, as servers do: build hg38 first, then hg19 with
--domains <hg38 out>/proteindomain.db, which adds the hg19 transcripts that the file lacks.

Rules, checked to reproduce the hg38 gene db of 2025-10-19 from the same sources:
- RefSeq rows: every row of ncbiRefSeqCurated; default isoform: the MANE Select transcript of the
  gene, else a member of ncbiRefSeqSelect (RefSeq Select); description: ncbiRefSeqLink.product
- GENCODE rows: every row of the comprehensive and pseudogene sets; default isoform and
  gene2canonicalisoform: the Ensembl_canonical tag, which mouse M25 does not have, so mm10 has no
  default GENCODE isoform; description: kgXref
- MANE (human): the NCBI summary of the current MANE release names each MANE Select and MANE Plus
  Clinical transcript as a versioned RefSeq/Ensembl pair. A gene model that is one of them, at
  that exact version, says so in its "mane" attribute. MANE is made on GRCh38: hg19 marks the
  transcripts whose version it still shares, and is not required to have them all
- retired RefSeq versions (hg38): from ncbiRefSeqHistorical, in a retiredgenes table. The server
  opens one only when asked for it by its exact accession
- protein domains (human): the features that NCBI's Conserved Domain Database annotates on a RefSeq
  protein, from the weekly RefSeq protein files, the only ones that carry them. Keyed by the
  transcript accession with its version. The positions are residue numbers, so a transcript gets
  the domains of its own protein only, matched on the protein accession and version that
  ncbiRefSeqLink gives it: the weekly files may be ahead of the annotation release that UCSC
  serves, and a transcript whose protein they have since revised gets no domains. A GENCODE
  transcript gets the domains of a RefSeq transcript with the same coding exons, which is the same
  protein. Features that NCBI copies from UniProt (modified residues, interaction regions) are
  left out
- genealias: the ENSG accession of a gene, and for human the HGNC alias and previous symbols
- codingGenes: genes with a coding RefSeq transcript; gene2coord: span of the default RefSeq
  isoform, or the widest one of a gene that has no default

Network: hgdownload.soe.ucsc.edu, ftp.ncbi.nlm.nih.gov (MANE, proteins), storage.googleapis.com
(HGNC), proteinpaint.stjude.org.
"""

import argparse, csv, datetime, gzip, json, os, re, sqlite3, sys, urllib.request
from email.utils import parsedate_to_datetime

UCSC = 'https://hgdownload.soe.ucsc.edu/goldenPath/{}/database/'
HGNC = 'https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt'
MANE = 'https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/'
PROTEINS = 'https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/mRNA_Prot/'
# the source database of a CDD model, as the last part of a feature note names it -> key of a domain row
DOMAIN_SOURCES = {'pfam': 'Pfam', 'smart': 'SMART', 'cd': 'Curated_at_NCBI', 'COG': 'COG', 'TIGR': 'TIGRFAM', 'PRK': 'PRK'}
OLD_DB = 'https://proteinpaint.stjude.org/ppSupport/genes.{}.db'
HUMAN = ['TP53', 'KRAS', 'EZH2', 'MYC', 'BRAF', 'EGFR', 'RUNX1', 'SETD2', 'KMT2D', 'BCR', 'ABL1']
"""what differs between genomes:
gencode    how a GENCODE release is spelled in the names of the UCSC tables, with the release number as the group
human      HGNC aliases, MANE and protein domains apply
mane       fewest MANE Select transcripts the build must find, per transcript set; 0 where MANE is not made on this assembly
retired    UCSC has ncbiRefSeqHistorical for this assembly
landmarks  genes that every build must have a default RefSeq isoform for, and a default GENCODE one where GENCODE marks any
refseq, ensembl   fewest rows of each transcript set
"""
GENOMES = {
	'hg38': dict(gencode=r'V(\d+)', human=True, mane=18000, retired=True, landmarks=HUMAN, refseq=90000, ensembl=300000),
	'hg19': dict(gencode=r'V(\d+)lift37', human=True, mane=0, retired=False, landmarks=HUMAN, refseq=80000, ensembl=200000),
	'mm10': dict(gencode=r'VM(\d+)', human=False, mane=0, retired=False, landmarks=['Trp53', 'Kras', 'Ezh2', 'Myc', 'Braf', 'Egfr', 'Runx1', 'Setd2', 'Kmt2d', 'Bcr', 'Abl1'], refseq=35000, ensembl=100000),
}

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


def latest_gencode(genome):
	with urllib.request.urlopen(UCSC.format(genome)) as r:
		versions = re.findall(r'wgEncodeGencodeComp' + GENOMES[genome]['gencode'] + r'\.txt\.gz', r.read().decode())
	if not versions:
		raise SystemExit('no wgEncodeGencodeComp table of this genome found at UCSC')
	return max(map(int, versions))


def clean(accession):
	"""ENST00000456328.2_1 (lifted to GRCh37) and ENST00000381192.10_PAR_Y (the chrY copy of a PAR transcript)
	are the accession ENST00000456328.2 and ENST00000381192.10, as RefSeq names a transcript wherever it lies"""
	return re.sub(r'(\.\d+)_\w+$', r'\1', accession)


def mane_release(cache):
	"""(release, {accession with version: 'Select' | 'Plus Clinical'}), for the RefSeq and the Ensembl transcript of each pair"""
	with urllib.request.urlopen(MANE) as r:
		found = re.search(r'MANE\.GRCh38\.v([\d.]+)\.summary\.txt\.gz', r.read().decode())
	if not found:
		raise SystemExit('no MANE summary file found at ' + MANE)
	status = {}
	for l in rows(fetch(MANE + found.group(0), cache)):
		if l[9].startswith('MANE '):
			status[l[5]] = status[l[7]] = l[9][5:]
	return found.group(1), status


def protein_features(path):
	"""a GenPept flat file -> (protein accession with version, protein length, [domain row]) of each protein, for the
	Region and Site features that carry a CDD cross reference"""

	def flush():
		if feature and 'CDD' in feature['q'] and feature['key'] in ('Region', 'Site'):
			q = feature['q']
			note, _, source = q.get('note', '').rpartition('; ')
			if feature['key'] == 'Site' or not re.fullmatch(r'[A-Za-z]+\d+', source):
				note, source = q.get('note', ''), ''
			row = {'name': q.get('region_name') or q.get('site_type') or 'other', 'description': note}
			# a domain names the model it is a hit of; a motif or site within one has no link of its own
			if source:
				row['CDD'] = q['CDD']
				key = DOMAIN_SOURCES.get(re.sub(r'\d+$', '', source))
				if key:
					row[key] = source
			# a site is one or more residues or ranges, each its own row
			for part in re.findall(r'(\d+)(?:\.\.>?(\d+))?', feature['loc']):
				record[2].append({**row, 'start': int(part[0]), 'stop': int(part[1] or part[0])})

	record = feature = name = None
	with gzip.open(path, 'rt', encoding='utf-8', errors='replace') as fh:
		for line in fh:
			if line.startswith('LOCUS'):
				record, feature = [None, int(line.split()[2]), []], None
			elif line.startswith('VERSION') and record:
				record[0] = line.split()[1]
			elif line.startswith('     ') and line[5] != ' ':
				flush()
				key, loc = line.split()[:2]
				feature, name = {'key': key, 'loc': loc, 'q': {}}, None
			elif feature and line.startswith(' ' * 21):
				text = line.strip()
				found = re.match(r'/(\w+)="?(.*)', text)
				if found:
					name, text = found.group(1), found.group(2)
					if name == 'db_xref' and text.startswith('CDD:'):
						name, text = 'CDD', text[4:]
					feature['q'][name] = text.rstrip('"')
				elif name == 'note' or name == 'region_name':  # a value continued on the next line
					feature['q'][name] = (feature['q'][name] + ' ' + text).rstrip('"')
			elif line.startswith(('ORIGIN', '//')) and record:
				flush()
				if record[0] and record[2]:
					yield tuple(record)
				record = feature = None


def build_domains(args, refseq, gencode, protein, today):
	"""writes the protein domain db, or adds to it the transcripts it lacks; returns a line for the report"""
	with urllib.request.urlopen(PROTEINS) as r:
		files = sorted(set(re.findall(r'human\.\d+\.protein\.gpff\.gz', r.read().decode())))
	if not files:
		raise SystemExit('no protein files found at ' + PROTEINS)
	paths = [fetch(PROTEINS + name, args.cache) for name in files]
	released = source_dates[files[0]]
	# the protein of a transcript: what its coding exons spell, so the same exons are the same protein
	cds = lambda chrom, m: (chrom, m['strand'], json.dumps(m['coding']))
	by_protein = {}
	for path in paths:
		for accession, aa, features in protein_features(path):
			by_protein[accession] = aa, [json.dumps(f, separators=(',', ':')) for f in features]
	domains = {}
	for chrom, _, _, m in refseq:
		aa, features = by_protein.get(protein.get(m['isoform']), (None, None))
		# also the protein length of the gene model, which differs for a model on an alt contig
		if 'cdslen' in m and aa == m['cdslen'] // 3 - 1 and is_primary(chrom):
			domains[m['isoform']] = features
	by_cds = {}
	for chrom, _, _, m in refseq:
		if m['isoform'] in domains and is_primary(chrom):
			by_cds.setdefault(cds(chrom, m), m['isoform'])
	copied = {}
	for chrom, _, _, m, _ in gencode:
		if 'cdslen' in m and is_primary(chrom) and cds(chrom, m) in by_cds:
			copied[m['isoform']] = domains[by_cds[cds(chrom, m)]]

	dbfile = args.domains or os.path.join(args.out, 'proteindomain.db')
	if not args.domains and os.path.exists(dbfile):
		os.remove(dbfile)
	db = sqlite3.connect(dbfile)
	db.executescript('CREATE TABLE IF NOT EXISTS domain (isoform varchar(30) not null, data json not null); CREATE TABLE IF NOT EXISTS buildDate (date char(100), sources text);')
	# a transcript of another genome's build is the same accession and version, so the same protein
	held = {i for (i,) in db.execute('select distinct isoform from domain')}
	db.executemany('insert into domain values (?,?)', [(i, d) for src in (domains, copied) for i, lst in src.items() if i not in held for d in lst])
	db.execute('insert into buildDate values (?,?)', (today, f'{args.genome}: NCBI CDD features of the RefSeq proteins of {released}'))
	db.execute('CREATE INDEX IF NOT EXISTS domain_isoform on domain(isoform collate nocase)')
	db.commit()
	assert len(domains) > 40000, f'only {len(domains)} RefSeq transcripts have domains'
	for gene in GENOMES[args.genome]['landmarks']:
		for src in (domains, copied):
			assert any(m['name'] == gene and m['isoform'] in src for g in (refseq, gencode) for m in (x[3] for x in g)), f'{gene} has no domains on any transcript'
	return f'Protein domains: NCBI RefSeq proteins of {released}; {len(domains)} RefSeq and {len(copied)} GENCODE transcripts.'


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
	genome, conf = args.genome, GENOMES[args.genome]
	v = args.gencode or latest_gencode(genome)
	gv = conf['gencode'].replace(r'(\d+)', str(v))  # e.g. V50, V45lift37, VM25
	release = gv[1:].replace('lift37', '')  # e.g. 50, 45, M25
	# in a dir of their own: the tables are named the same for every genome
	os.makedirs(os.path.join(args.cache, genome), exist_ok=True)
	f = lambda name: fetch(UCSC.format(genome) + name + '.txt.gz', os.path.join(args.cache, genome))
	old_db = args.old_db or fetch(OLD_DB.format(genome), args.cache)

	# ---------- RefSeq
	product, protein = {}, {}  # description and protein accession of each transcript
	for l in rows(f('ncbiRefSeqLink')):
		product[l[0]], protein[l[0]] = l[3], l[5]
	mane_version, mane = mane_release(args.cache) if conf['human'] else (None, {})
	select = {l[1] for l in rows(f('ncbiRefSeqSelect'))}
	# a gene's MANE Select transcript is its default; RefSeq Select decides only for a gene with none
	mane_genes = {l[12] for l in rows(f('ncbiRefSeqCurated')) if mane.get(l[1]) == 'Select'}
	select = {l[1] for l in rows(f('ncbiRefSeqCurated')) if mane.get(l[1]) == 'Select' or (l[1] in select and l[12] not in mane_genes)}
	refseq = [gene_model(l, l[12], product.get(l[1])) for l in rows(f('ncbiRefSeqCurated'))]
	for _, _, _, m in refseq:
		if m['isoform'] in mane:
			m['mane'] = mane[m['isoform']]

	# ---------- GENCODE
	attrs = {clean(l[4]): l for l in rows(f(f'wgEncodeGencodeAttrs{gv}'))}  # transcript accession -> row
	canonical = {clean(l[0]) for l in rows(f(f'wgEncodeGencodeTag{gv}')) if l[1] == 'Ensembl_canonical'}
	# keyed with and without the version: kgXref is not always of the GENCODE release being built
	kgdesc = {}
	for l in rows(f('kgXref')):
		kgdesc[clean(l[0])] = l[7]
		kgdesc.setdefault(base(clean(l[0])), l[7])
	gencode = []
	for l in [*rows(f(f'wgEncodeGencodeComp{gv}')), *rows(f(f'wgEncodeGencodePseudoGene{gv}'))]:
		l[1] = clean(l[1])
		a = attrs.get(l[1])
		desc, category = kgdesc.get(l[1]) or kgdesc.get(base(l[1])), a[12] if a else None
		gencode.append(gene_model(l, l[12], desc, category) + (l[1] in canonical,))
		if l[1] in mane:
			gencode[-1][3]['mane'] = mane[l[1]]

	write_track(list(refseq), os.path.join(args.out, f'refGene.{genome}.gz'))
	write_track([g[:4] for g in gencode], os.path.join(args.out, f'gencode.v{release}.{genome}.gz'))

	# ---------- gene db
	dbfile = os.path.join(args.out, f'genes.{genome}.db')
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

	genes = [db_row(*g, g[3]['isoform'] in select, 'RefSeq') for g in refseq]
	genes += [db_row(*g[:4], g[4], 'GENCODE') for g in gencode]
	# the rows of the previous db that no source table has: hand-made models named other than by an accession,
	# such as the chrM polycistronic transcript
	old = sqlite3.connect(f'file:{old_db}?mode=ro', uri=True)
	genes += [r for r in old.execute('select name,isoform,isdefault,genemodel from genes') if r[1] and not re.match(r'[NXY][MRP]_|ENS', r[1])]
	db.executemany('insert into genes values (?,?,?,?)', genes)

	if conf['retired']:
		current = {m['isoform'] for _, _, _, m in refseq}
		retired = [gene_model(l, l[12]) for l in rows(f('ncbiRefSeqHistorical')) if l[1] not in current]
		retired = [db_row(c, a, b, {**m, 'retired': True}, 0, 'RefSeq') for c, a, b, m in retired]
		db.executescript('CREATE TABLE retiredgenes (name character varying(50) collate nocase, isoform character varying(50) collate nocase, genemodel json);')
		db.executemany('insert into retiredgenes values (?,?,?)', [(r[0], r[1], r[3]) for r in retired])
		db.executescript('CREATE INDEX retiredisoform on retiredgenes(isoform collate nocase);')

	alias = {(base(clean(a[0])), a[1]) for a in attrs.values()}  # ENSG -> symbol
	if conf['human']:
		with open(fetch(HGNC, args.cache), encoding='utf-8') as fh:
			for r in csv.DictReader(fh, delimiter='\t'):
				for key in ('alias_symbol', 'prev_symbol'):
					alias.update((a, r['symbol']) for a in r[key].strip('"').split('|') if a)
				if r['ensembl_gene_id']:
					alias.add((r['ensembl_gene_id'], r['symbol']))
	db.executemany('insert into genealias values (?,?)', sorted(alias))

	# ENST unversioned, as the consumer queries the gdc api with it; refseq in both forms
	r2e = set()
	for l in rows(f(f'wgEncodeGencodeRefSeq{gv}')):
		enst = base(clean(l[0]))
		r2e.update({(enst, l[1]), (enst, base(l[1]))})
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
				g2c.setdefault(base(clean(a[0])), m['isoform'])
	db.executemany('insert into gene2canonicalisoform values (?,?)', g2c.items())

	coding = sorted({m['name'] for _, _, _, m in refseq if 'cdslen' in m})
	db.executemany('insert into codingGenes values (?)', [(n,) for n in coding])

	today = datetime.date.today().isoformat()
	sources = f"RefSeq: UCSC {genome} ncbiRefSeqCurated {source_dates['ncbiRefSeqCurated.txt.gz']}; GENCODE {'v' + release}{' lifted to GRCh37' if 'lift37' in gv else ''}"
	if conf['human']:
		sources += f"; MANE v{mane_version}; HGNC {source_dates['hgnc_complete_set.txt']}"
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

	selfcheck(db, conf, bool(canonical))
	if conf['human']:
		sources += '\n\n' + build_domains(args, refseq, gencode, protein, today)
	diff(old, db, args.out, today, sources, genome)
	print(f'built {dbfile}: {sources}', file=sys.stderr)


def selfcheck(db, conf, has_canonical):
	"""stop the build rather than publish a db that would break a server"""
	q = lambda sql, *p: db.execute(sql, p).fetchone()[0]
	sets = ('NM_%', 'ENS%') if has_canonical else ('NM_%',)  # the transcript sets that mark a default isoform
	assert q("select count(*) from genes where isoform like 'N%'") > conf['refseq'], 'too few RefSeq rows'
	assert q("select count(*) from genes where isoform like 'ENS%'") > conf['ensembl'], 'too few GENCODE rows'
	assert q("select count(*) from genes where (isoform like 'N%' or isoform like 'ENS%') and isoform not glob '*.[0-9]*'") == 0, 'accession without version'
	assert q("select count(*) from genes where isoform glob '*.[0-9]*_*'") == 0, 'accession with a suffix after its version'
	assert q('select count(*) from genealias') > (100000 if conf['human'] else 40000), 'too few aliases'
	if conf['retired']:
		assert q('select count(*) from retiredgenes') > 50000, 'too few retired versions'
		assert q('select count(*) from retiredgenes r join genes g on g.isoform=r.isoform') == 0, 'a version is both current and retired'
	assert q('select count(*) from codingGenes') > 19000, 'too few coding genes'
	if conf['mane']:
		for prefix in sets:
			n = q("""select count(*) from genes where isoform like ? and json_extract(genemodel,'$.mane')='Select'""", prefix)
			assert n > conf['mane'], f'only {n} {prefix[:-1]} transcripts are MANE Select: the MANE release does not match this annotation'
	assert q("""select count(*) from genes where isoform like 'N%' and json_extract(genemodel,'$.mane')='Select' and isdefault=0""") == 0, 'a MANE Select RefSeq transcript is not the default'
	for gene in conf['landmarks']:
		for prefix in sets:
			n = q('select count(*) from genes where name=? and isdefault=1 and isoform like ?', gene, prefix)
			assert n >= 1, f'{gene} has no default {prefix[:-1]} isoform'
			if conf['mane']:
				n = q("""select count(*) from genes where name=? and isdefault=1 and isoform like ? and json_extract(genemodel,'$.mane')='Select'""", gene, prefix)
				assert n >= 1, f'the default {prefix[:-1]} isoform of {gene} is not its MANE Select'
		assert q('select count(*) from gene2coord where name=?', gene) == 1, f'{gene} missing from gene2coord'
		if has_canonical:
			assert q('select count(*) from gene2canonicalisoform where gene=?', gene) == 1, f'{gene} has no canonical isoform'


def diff(old, new, out, today, sources, genome):
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
	dated = old.execute("select count(*) from sqlite_master where name='buildDate'").fetchone()[0]
	old_date = old.execute('select date from buildDate').fetchone()[0] if dated else 'at an unrecorded date'
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
		fh.write(f'# {genome} gene models {today}\n\n{sources}\n\nCompared with the gene db built {old_date}, primary chromosomes only.\n\n| Change | Count |\n|---|---:|\n')
		fh.writelines(f'| {labels[k]} | {len(v)} |\n' for k, v in lists.items())
		for k in ('protein_changed', 'symbol_renamed', 'default_changed'):
			fh.write(f'\n## {labels[k].lstrip(". ")} (first 50; all in diff.json)\n\n')
			fh.writelines(f'- {i}\n' for i in lists[k][:50])


if __name__ == '__main__':
	p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
	p.add_argument('--genome', default='hg38', choices=list(GENOMES))
	p.add_argument('--out', required=True, help='output dir')
	p.add_argument('--domains', help='existing protein domain db to add the transcripts of this genome to, default: write a new <out>/proteindomain.db')
	p.add_argument('--gencode', type=int, help='GENCODE release number, default: latest at UCSC for the genome')
	p.add_argument('--old-db', help='previous gene db of the genome, default: download the published one')
	p.add_argument('--cache', help='dir for downloads, default: <out>/downloads')
	args = p.parse_args()
	args.cache = args.cache or os.path.join(args.out, 'downloads')
	build(args)
