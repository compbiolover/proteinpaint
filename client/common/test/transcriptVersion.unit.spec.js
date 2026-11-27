import tape from 'tape'
import { isoformLabel, showsVersions, setShowVersions } from '../transcriptVersion.js'

/*
Tests whether a transcript accession is shown with its version

	isoformLabel()
	showsVersions(), setShowVersions()
*/

tape('\n', test => {
	test.comment(`-***- common/transcriptVersion unit -***-`)
	test.end()
})

tape('isoformLabel()', test => {
	const was = showsVersions()

	setShowVersions(false)
	test.equal(showsVersions(), false, 'versions are hidden once turned off')
	test.equal(isoformLabel('NM_004456.5'), 'NM_004456', 'an accession is shown without its version')
	test.equal(isoformLabel({ isoform: 'ENST00000320356.7' }), 'ENST00000320356', 'so is that of a gene model')
	test.equal(isoformLabel('NM_004456'), 'NM_004456', 'an accession with no version is shown as is')
	test.equal(
		isoformLabel({ isoform: 'NM_004456.1', retired: true }),
		'NM_004456.1',
		'a retired version always shows its version'
	)

	setShowVersions(true)
	test.equal(showsVersions(), true, 'versions are shown once turned on')
	test.equal(isoformLabel('NM_004456.5'), 'NM_004456.5', 'an accession is shown with its version')
	test.equal(isoformLabel({ isoform: 'ENST00000320356.7' }), 'ENST00000320356.7', 'so is that of a gene model')

	setShowVersions(was)
	test.end()
})
