import { stripAccessionVersion } from '#shared/helpers.js'

/*
Whether transcript accessions are shown with their version (NM_004456.5) or without (NM_004456).
Without by default: the version is what a lookup or a comparison of two releases needs, and noise
to a user who only wants to see the gene. A preference of this browser, as a user who turns it on
wants it on in every view.
*/

const key = 'sjpp-show-transcript-versions'

export function showsVersions() {
	try {
		return localStorage.getItem(key) == '1'
	} catch (e) {
		return false // storage is blocked: the default
	}
}

export function setShowVersions(on) {
	try {
		localStorage.setItem(key, on ? '1' : '0')
	} catch (e) {
		// storage is blocked: the choice lasts until the next read, which is none
	}
}

/* The accession of a gene model, or an accession, as shown to the user. A retired version always
shows its version, as that is what tells it from the current transcript. */
export function isoformLabel(gm) {
	const isoform = typeof gm == 'string' ? gm : gm.isoform
	return showsVersions() || gm.retired ? isoform : stripAccessionVersion(isoform)
}
