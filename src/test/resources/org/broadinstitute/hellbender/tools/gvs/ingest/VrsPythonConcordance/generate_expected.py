#!/usr/bin/env python3
"""Generate expected.tsv, the vrs-python IDs that VrsPythonConcordanceUnitTest compares GVS against.

Each case is a VCF record as GVS ingests it, written against the GRCh38 analysis set
(src/test/resources/large/Homo_sapiens_assembly38.20.21.fasta). vrs-python normalizes and identifies
every ALT against the RefSeq sequence the refget accession names, as any other VRS implementation
would. That is the point of the comparison: on chr21 the analysis set N-masks sequence that RefSeq
has, so GVS and vrs-python see different reference bases there (VS-2059).

Requires vrs-python 2.3.3, pinned below, with its extras, since the translator imports hgvs. It also
needs network access to NCBI unless --refseq-dir already holds the sequences:

    python3 -m venv vrsenv && vrsenv/bin/pip install 'ga4gh.vrs[extras]==2.3.3'
    vrsenv/bin/python generate_expected.py --refseq-dir /tmp/refseq > expected.tsv

If the download fails, for example with CERTIFICATE_VERIFY_FAILED, fetch the files with any other tool
into --refseq-dir as <accession>.fa, from
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=<accession>&rettype=fasta

No SeqRepo is needed. The script fetches NC_000020.11 and NC_000021.9 from NCBI E-utilities, checks
that each one's GA4GH digest equals the accession GVS uses (RefGetAccessionProvider), and serves
them to vrs-python through an in-memory data proxy.
"""

import argparse
import base64
import hashlib
import importlib.metadata
import sys
import urllib.request
from pathlib import Path

VRS_PYTHON_VERSION = "2.3.3"

# contig as in the GVS FASTA -> (RefSeq accession, refget accession). Must agree with RefGetAccessionProvider.
SEQUENCES = {
    "chr20": ("NC_000020.11", "SQ.-A1QmD_MatoqxvgVxBLZTONHz9-c7nQo"),
    "chr21": ("NC_000021.9", "SQ.5ZUqxCmDDgN4xTRbaSjN8LwgZironmB8"),
}

# Spanning deletion. VRS cannot represent it, so GVS writes an empty-string placeholder ID.
SPAN_DEL = "*"

# (case, contig, 1-based pos, REF, ALTs, note). REF is as the analysis set has it, so N inside a masked region.
# gvs_differs, when set, records why GVS deliberately does not match vrs-python for that ALT.
CASES = [
    ("snv", "chr20", 1000000, "C", ["T"], None),
    ("mnv", "chr20", 1005008, "GC", ["AT"], None),
    # Inserted TC follows anchor G, then GCTTG: no rotation of TC matches either flank, so nothing rolls.
    ("insertion_novel", "chr20", 1005007, "G", ["GTC"], None),
    # One more CA in the (CA)5 run at chr20:1045534-1045543: reference-derived, so a ReferenceLengthExpression.
    ("insertion_repeat", "chr20", 1045533, "T", ["TCA"], None),
    # Deletes GCT from ...GGG GCT TGT...; neither flank continues the deleted sequence, so nothing rolls.
    ("deletion_unique", "chr20", 1005007, "GGCT", ["G"], None),
    # Deletes one CCA from (CCA)4 at chr20:1012545-1012556, written left-aligned...
    ("deletion_repeat", "chr20", 1012544, "GCCA", ["G"], None),
    # ...and the same deletion written at the right end of the repeat. It must get the same IDs.
    ("deletion_repeat_right_shifted", "chr20", 1012553, "ACCA", ["A"], None),
    # Delins: both sides non-empty after trimming, so VRS does not roll it. CCA -> T inside the (CCA)4 above...
    ("delins_shorter_alt", "chr20", 1012544, "GCCA", ["GT"], None),
    # ...and CA -> GGGG inside the (CA)5 run of insertion_repeat. Both stay a LiteralSequenceExpression at the trimmed interval.
    ("delins_longer_alt", "chr20", 1045533, "TCA", ["TGGGG"], None),
    # TA(6)T at chr20:1005240-1005247: two deletions, an insertion and an SNV of the anchor, plus *.
    ("multiallelic", "chr20", 1005240, "TAA", ["T", "TA", SPAN_DEL, "TAAA", "GAA"], None),
    # chr21:11000000 is N in the analysis set. An SNV is unaffected, as long as ALT is not the RefSeq base.
    ("masked_snv", "chr21", 11000000, "N", ["C"], None),
    # Here ALT is the RefSeq base. vrs-python sees no change and returns a reference allele; GVS reads N,
    # sees a substitution, and returns a LiteralSequenceExpression. GVS normalizes against its own reference.
    ("masked_snv_matches_refseq", "chr21", 11000000, "N", ["A"],
     "analysis set N-masks chr21:11000000; RefSeq has A there, so vrs-python sees a reference allele"),
]

COLUMNS = ["case", "contig", "pos", "ref", "alts", "alt_index", "alt", "allele_id", "location_id",
           "start", "end", "state_type", "state_sequence", "state_length", "state_repeat_subunit_length",
           "gvs_differs"]


def sha512t24u(seq):
    return base64.urlsafe_b64encode(hashlib.sha512(seq.encode("ascii")).digest()[:24]).decode("ascii")


def read_fasta(path):
    lines = Path(path).read_text().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"{path} is not a FASTA file")
    return "".join(line.strip() for line in lines[1:]).upper()


def load_refseq(refseq_dir):
    refseq_dir.mkdir(parents=True, exist_ok=True)
    sequences = {}
    for contig, (accession, refget) in SEQUENCES.items():
        path = refseq_dir / f"{accession}.fa"
        if not path.exists():
            url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta"
            print(f"downloading {accession}", file=sys.stderr)
            urllib.request.urlretrieve(url, path)
        seq = read_fasta(path)
        digest = "SQ." + sha512t24u(seq)
        if digest != refget:
            raise ValueError(f"{path}: digest {digest} is not {refget}")
        sequences[refget] = (contig, accession, seq)
    return sequences


def make_data_proxy(sequences):
    from ga4gh.vrs.dataproxy import _DataProxy

    def refget_of(identifier):
        for refget, (contig, accession, _) in sequences.items():
            if identifier in (f"ga4gh:{refget}", refget, f"refseq:{accession}",
                              f"GRCh38:{contig}", f"GRCh38:{contig.removeprefix('chr')}"):
                return refget
        raise KeyError(identifier)

    class InMemoryDataProxy(_DataProxy):
        def get_sequence(self, identifier, start=None, end=None):
            return sequences[refget_of(identifier)][2][start:end]

        def get_metadata(self, identifier):
            refget = refget_of(identifier)
            contig, accession, seq = sequences[refget]
            return {"length": len(seq), "alphabet": "ACGT",
                    "aliases": [f"ga4gh:{refget}", f"refseq:{accession}",
                                f"GRCh38:{contig}", f"GRCh38:{contig.removeprefix('chr')}"]}

    return InMemoryDataProxy()


def vrs_allele(data_proxy, refget, pos, ref, alt):
    """Normalize and identify one ALT the way GVS sees it: interval [pos - 1, pos - 1 + len(REF)), state ALT."""
    from ga4gh.core import ga4gh_identify
    from ga4gh.vrs import models
    from ga4gh.vrs.normalize import normalize

    allele = models.Allele(
        location=models.SequenceLocation(
            sequenceReference=models.SequenceReference(refgetAccession=refget),
            start=pos - 1, end=pos - 1 + len(ref)),
        state=models.LiteralSequenceExpression(sequence=alt))
    allele = normalize(allele, data_proxy)
    allele.location.id = ga4gh_identify(allele.location)
    allele.id = ga4gh_identify(allele)
    return allele


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--refseq-dir", type=Path, required=True,
                        help="directory holding NC_000020.11.fa and NC_000021.9.fa; downloaded there if absent")
    args = parser.parse_args()

    installed = importlib.metadata.version("ga4gh.vrs")
    if installed != VRS_PYTHON_VERSION:
        sys.exit(f"vrs-python {VRS_PYTHON_VERSION} is required, found {installed}")
    from ga4gh.vrs.extras.translator import AlleleTranslator

    sequences = load_refseq(args.refseq_dir)
    refget_by_contig = {contig: refget for contig, (_, refget) in SEQUENCES.items()}
    data_proxy = make_data_proxy(sequences)
    translator = AlleleTranslator(data_proxy=data_proxy)

    out = sys.stdout
    out.write(f"# Generated by generate_expected.py with vrs-python (ga4gh.vrs) {VRS_PYTHON_VERSION} against RefSeq "
              f"{', '.join(a for a, _ in SEQUENCES.values())}. Do not edit by hand.\n")
    out.write("\t".join(COLUMNS) + "\n")
    for case, contig, pos, ref, alts, gvs_differs in CASES:
        refget = refget_by_contig[contig]
        refseq_ref = sequences[refget][2][pos - 1:pos - 1 + len(ref)]
        if "N" not in ref and ref != refseq_ref:
            raise ValueError(f"{case}: REF {ref} is not the RefSeq sequence {refseq_ref}")
        for alt_index, alt in enumerate(alts):
            row = {"case": case, "contig": contig, "pos": pos, "ref": ref, "alts": ",".join(alts),
                   "alt_index": alt_index, "alt": alt, "gvs_differs": gvs_differs or ""}
            if alt == SPAN_DEL:
                rows = {c: "" for c in COLUMNS}
                rows.update(row)
                row = rows
            else:
                allele = vrs_allele(data_proxy, refget, pos, ref, alt)
                if "N" not in ref:
                    # The translator is how most callers reach vrs-python; it must agree with the direct path.
                    translated = translator.translate_from(f"{contig.removeprefix('chr')}-{pos}-{ref}-{alt}", "gnomad")
                    if translated.id != allele.id:
                        raise ValueError(f"{case} {alt}: translator gives {translated.id}, normalize gives {allele.id}")
                state = allele.state
                row.update({
                    "allele_id": allele.id,
                    "location_id": allele.location.id,
                    "start": allele.location.start,
                    "end": allele.location.end,
                    "state_type": state.type,
                    "state_sequence": state.sequence.root if state.type == "LiteralSequenceExpression" else "",
                    "state_length": state.length if state.type == "ReferenceLengthExpression" else "",
                    "state_repeat_subunit_length":
                        state.repeatSubunitLength if state.type == "ReferenceLengthExpression" else "",
                })
            out.write("\t".join(str(row[c]) for c in COLUMNS) + "\n")


if __name__ == "__main__":
    main()
