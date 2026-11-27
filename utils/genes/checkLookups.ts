/*
Runs the server's gene lookups against two real gene dbs and prints, for well known transcripts
asked for without a version, at the current version and at a retired one, which model each db
opens. Exits 1 when an accession does not resolve in the new db, or resolves to another gene in
either. For a human genome: the transcripts are those of human genes.
A protein length that differs between the two dbs is printed for review, not failed: a db that
holds retired versions opens the retired model, and a new release may change a current one.

	cd server && npx tsx --conditions=sjpp/dev ../utils/genes/checkLookups.ts <old genes db> <new genes db>
*/
import Database from 'better-sqlite3'
import { initGeneDbLookups } from '../../server/src/genedbLookups.ts'
import { getResult } from '../../server/src/gene.js'

function load(file: string) {
	const db = new Database(file, { readonly: true, fileMustExist: true })
	const genedb: any = {
		db,
		getjsonbyname: db.prepare('select isdefault,genemodel from genes where name=?'),
		getjsonbyisoform: db.prepare('select isdefault,genemodel from genes where isoform=?')
	}
	const tables = db.prepare("select name from sqlite_master where type='table'").all() as { name: string }[]
	initGeneDbLookups(genedb, new Set(tables.map(t => t.name)))
	return { genedb, genomicNameRegexp: /[^a-zA-Z0-9.:_-]/ }
}

// [gene, accessions to ask for]; a retired version opens its own model in a db that holds retired versions
const cases: [string, string[]][] = [
	[
		'EZH2',
		['NM_004456', 'NM_004456.1', 'NM_004456.5', 'ENST00000320356', 'ENST00000320356.7', 'ENSG00000106462', 'KMT6']
	],
	[
		'TP53',
		['NM_000546', 'NM_000546.5', 'NM_000546.6', 'ENST00000269305', 'enst00000269305.4', 'ENSG00000141510', 'p53']
	],
	['MYC', ['NM_002467', 'NM_002467.2', 'NM_002467.6']],
	['RUNX1', ['NM_001754', 'NM_001754.1', 'NM_001754.5']],
	['SETD2', ['NM_014159', 'NM_014159.5', 'NM_014159.7']],
	['KMT2D', ['NM_003482', 'NM_003482.2', 'NM_003482.4']],
	['KRAS', ['NM_004985', 'NM_004985.5', 'NM_033360', 'ENST00000311936']],
	['BRAF', ['NM_004333', 'NM_004333.4', 'ENST00000646891']],
	['EGFR', ['NM_005228', 'NM_005228.3']]
]

/** the model a block would open for this input: the isoform asked for, else the gene's default */
function opened(g: any, input: string) {
	const r: any = getResult(g, { input, deep: true })
	const primary = r.gmlst.filter((m: any) => !m.chr.includes('_'))
	const asked = (r.found_isoform || input).toUpperCase()
	return primary.find((m: any) => m.isoform.toUpperCase() == asked) || primary.find((m: any) => m.isdefault)
}
const show = (m: any) =>
	m ? `${m.isoform} ${m.cdslen ? Math.floor(m.cdslen / 3) - 1 + ' aa' : 'noncoding'}` : 'NOT FOUND'

const [oldFile, newFile] = process.argv.slice(2)
if (!oldFile || !newFile) throw 'usage: checkLookups.ts <old genes db> <new genes db>'
const oldDb = load(oldFile)
const newDb = load(newFile)
let failed = 0
for (const [gene, inputs] of cases) {
	for (const input of inputs) {
		const a = opened(oldDb, input)
		const b = opened(newDb, input)
		// the new db must find it; the old one may not have known the name at all, e.g. an ENSG before aliases held them
		const ok = b?.name == gene && (!a || a.name == gene)
		if (!ok) failed++
		const note = ok && a && a.cdslen != b.cdslen ? '  <- protein length differs, review' : ''
		console.log(
			`${ok ? 'ok  ' : 'FAIL'} ${gene.padEnd(6)} ${input.padEnd(19)} old: ${show(a).padEnd(26)} new: ${show(b)}${note}`
		)
	}
}
for (const g of [oldDb, newDb]) {
	if (g.genedb.resolveIsoform('NM_999999.1') || getResult(g, { input: 'NOTAGENE', deep: true }).gmlst.length) {
		failed++
		console.log('FAIL a made-up name resolved')
	}
}
console.log(failed ? `${failed} FAILED` : 'all resolved')
process.exit(failed ? 1 : 0)
