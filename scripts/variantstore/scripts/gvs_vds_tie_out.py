import hail as hl
import pytest

mt = hl.read_matrix_table('joined.mt')
rs = mt.rows()
cs = mt.cols()
es = mt.entries()


"""
Schema of `mt`:
----------------------------------------
Global fields:
    'tranche_data': array<struct {
        model: str,
        truth_sensitivity: float64,
        min_vqslod: float64,
        filter_name: str
    }>
    'truth_sensitivity_snp_threshold': float64
    'truth_sensitivity_indel_threshold': float64
    'snp_vqslod_threshold': float64
    'indel_vqslod_threshold': float64
----------------------------------------
Column fields:
    's': str
    'vcf_col': struct {
        s: str
    }
    'vds_col': struct {
        s: str
    }
----------------------------------------
Row fields:
    'locus': locus<GRCh38>
    'vcf_row': struct {
        locus: locus<GRCh38>,
        alleles: array<str>,
        rsid: str,
        qual: float64,
        filters: set<str>,
        info: struct {
            AC: array<int32>,
            AF: array<float64>,
            AN: int32,
            AS_QUALapprox: str,
            AS_VQSLOD: array<str>,
            AS_YNG: array<str>,
            QUALapprox: int32
        }
    }
    'vds_row': struct {
        locus: locus<GRCh38>,
        alleles: array<str>,
        filters: set<str>,
        as_vqsr: dict<str, struct {
            model: str,
            vqslod: float64,
            yng_status: str
        }>
    }
----------------------------------------
Entry fields:
    'vcf_entry': struct {
        AD: array<int32>,
        FT: str,
        GQ: int32,
        GT: call,
        RGQ: int32
    }
    'vds_entry': struct {
        LGT: call,
        GQ: int32,
        RGQ: int32,
        LAD: array<int32>,
        LA: array<int32>,
        FT: bool
    }
----------------------------------------
Column key: ['s']
Row key: ['locus']
----------------------------------------
"""


def updated_local_alleles(mt=mt):
    """Creates field vds_entry.LA_vcf reindexing the VDS LA field to the VCF allele order"""
    mt2 = mt.annotate_rows(vcf_allele_dict=hl.dict(hl.enumerate(mt.vcf_row.alleles, index_first=False)))
    mt2 = mt2.annotate_rows(vds_allele_mapping=mt2.vds_row.alleles.map(lambda a: mt2.vcf_allele_dict[a]))
    mt2 = mt2.annotate_entries(
        vds_entry=mt2.vds_entry.annotate(LA_vcf=mt2.vds_entry.LA.map(lambda i: mt2.vds_allele_mapping[i])))
    return mt2


def gq0_filter(mt=mt):
    """Sets VDS LGTs to missing at reference genotypes when GQ is 0"""
    mt2 = mt.annotate_entries(vds_entry=mt.vds_entry.annotate(
        LGT=hl.or_missing(mt.vds_entry.LGT.is_non_ref() | (mt.vds_entry.GQ > 0), mt.vds_entry.LGT)))
    return mt2


def ft_filter(mt=mt):
    """Sets VDS LGTs to missing at non-ref genotypes when FT is False"""
    mt2 = mt.annotate_entries(
        vds_entry=mt.vds_entry.annotate(LGT=hl.or_missing(hl.coalesce(mt.vds_entry.FT, True), mt.vds_entry.LGT)))
    return mt2


def find_violations(ht, predicate):
    """Checks that `predicate` is true for all rows of `ht`.
    Parameters
    ----------
    ht
        Hail table.
    predicate
        Boolean expression in the row index of the `ht`.
    """
    ht = ht.annotate(CHECK=predicate)
    violations = ht.aggregate(hl.agg.count_where(~ht.CHECK))
    if violations > 0:
        ht = ht.filter(~ht.CHECK)
        ht.show()
        assert False, f'{violations} violations'


def test_locus_equivalence():
    find_violations(rs, hl.is_defined(rs.vcf_row) & hl.is_defined(rs.vds_row))


def test_sample_equivalence():
    find_violations(cs, hl.is_defined(cs.vcf_col) & hl.is_defined(cs.vds_col))


def test_allele_equivalence():
    find_violations(rs, hl.set(rs.vcf_row.alleles) == hl.set(rs.vds_row.alleles))


def test_filter_equivalence():
    find_violations(rs, rs.vcf_row.filters == rs.vds_row.filters)


def test_gt_equivalence():
    es = updated_local_alleles().entries()
    find_violations(es, hl.vds.lgt_to_gt(es.vds_entry.LGT, hl.coalesce(es.vds_entry.LA_vcf, [0])) == es.vcf_entry.GT)


def test_gq_equivalence():
    find_violations(es, es.vds_entry.GQ == es.vcf_entry.GQ)


def test_ad_equivalence():
    es = updated_local_alleles().entries()

    def lad_to_ad(lad, la, n_alleles):
        return hl.rbind(hl.dict(hl.zip(la, lad)), lambda d: hl.range(n_alleles).map(lambda i: d.get(i, 0)))

    es = es.annotate(
        vds_AD=lad_to_ad(es.vds_entry.LAD, hl.coalesce(es.vds_entry.LA_vcf, [0]), hl.len(es.vcf_row.alleles)))
    find_violations(es, es.vds_AD == es.vcf_entry.AD)


def test_ft_equivalence():
    es = mt.filter_rows(hl.len(mt.vcf_row.filters) == 0).entries()
    es2 = es.filter(es.vcf_entry.GT.is_non_ref())
    ft1 = es2.vds_entry.FT
    ft2 = es2.vcf_entry.FT
    find_violations(es2, hl.case()
                    .when(hl.is_missing(ft1) != hl.is_missing(ft2), False)
                    .when((ft1 == True) & (ft2 == 'PASS'), True)
                    .when((ft1 == False) & (ft2 != 'PASS'), True)
                    .default(False))


def test_call_stats_equivalence_an():
    mt = updated_local_alleles(ft_filter(gq0_filter()))
    mt = mt.annotate_rows(
        vds_call_stats=hl.agg.call_stats(hl.vds.lgt_to_gt(mt.vds_entry.LGT, hl.coalesce(mt.vds_entry.LA_vcf, [0])),
                                         hl.len(mt.vds_row.alleles)))
    r = mt.rows()
    find_violations(r, (r.vds_call_stats.AN == r.vcf_row.info.AN))


def test_call_stats_equivalence_ac():
    mt = updated_local_alleles(ft_filter(gq0_filter()))
    mt = mt.annotate_rows(
        vds_call_stats=hl.agg.call_stats(hl.vds.lgt_to_gt(mt.vds_entry.LGT, hl.coalesce(mt.vds_entry.LA_vcf, [0])),
                                         hl.len(mt.vds_row.alleles)))
    r = mt.rows()
    find_violations(r, (r.vds_call_stats.AC[1:] == r.vcf_row.info.AC))


def test_call_stats_equivalence_af():
    mt = updated_local_alleles(ft_filter(gq0_filter()))
    mt = mt.annotate_rows(
        vds_call_stats=hl.agg.call_stats(hl.vds.lgt_to_gt(mt.vds_entry.LGT, hl.coalesce(mt.vds_entry.LA_vcf, [0])),
                                         hl.len(mt.vds_row.alleles)))
    r = mt.rows()
    find_violations(r, hl.all(lambda x: x, hl._zip_func(r.vds_call_stats.AF[1:], r.vcf_row.info.AF,
                                                        f=lambda t1, t2: hl.approx_equal(t1, t2, tolerance=0.01))))
