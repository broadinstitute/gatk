from scorevariants.readers import (
    determine_polymorphic_type,
    ReferenceTensorReader,
    TensorReader,
)
from scorevariants.utilities import (
    variant_is_snp,
)
from scorevariants.encoders import (
    AnnotationEncoder,
    ReferenceEncoder,
)

from scorevariants.tests.utilities import MockVariantRecord
import numpy.testing as npt


def test_determine_polymoprhic_type():
    hom_snp = MockVariantRecord(contig="chr1", pos=14396, ref="A", alts=("C",))
    assert determine_polymorphic_type(hom_snp) == "SNP"

    het_mnp = MockVariantRecord(
        contig="chr1", pos=14396, ref="TA", alts=("CG", "TG")
    )
    assert determine_polymorphic_type(het_mnp) == "MNP"

    het_indel = MockVariantRecord(
        contig="chr1", pos=14396, ref="TA", alts=("C", "T")
    )
    assert determine_polymorphic_type(het_indel) == "INDEL"

    mixed_snp_indel = MockVariantRecord(
        contig="chr1", pos=14396, ref="T", alts=("C", "TG")
    )
    assert determine_polymorphic_type(mixed_snp_indel)


def test_variant_is_snp():

    hom_snp = MockVariantRecord(contig="chr1", pos=14396, ref="A", alts=("C",))
    assert variant_is_snp(hom_snp)

    het_snp = MockVariantRecord(
        contig="chr1", pos=14396, ref="A", alts=("C", "*")
    )

    assert variant_is_snp(het_snp)

    hom_insert = MockVariantRecord(
        contig="chr1", pos=14396, ref="A", alts=("CAGCT",)
    )
    assert not variant_is_snp(hom_insert)

    het_insert = MockVariantRecord(
        contig="chr1", pos=14396, ref="A", alts=("C", "GGT")
    )
    assert not variant_is_snp(het_insert)

    hom_del = MockVariantRecord(
        contig="chr1", pos=14396, ref="AGCCT", alts=("A")
    )
    assert not variant_is_snp(hom_del)

    het_del = MockVariantRecord(
        contig="chr1", pos=14396, ref="AGCCT", alts=("A", ".")
    )
    assert not variant_is_snp(het_del)


seq_dict = {
    "contig01": "ACGTTCG",
    "contig02": "TGCATGC",
}


def test_reference_tensor_reader():
    reference_encoder = ReferenceEncoder(window=4, handle_start="pad")
    annotation_encoder = AnnotationEncoder(annotation_list=["MQ", "DP"])

    reference_reader = ReferenceTensorReader(
        variant_file=None,
        reference_dict=seq_dict,
        reference_encoder=reference_encoder,
        annotation_encoder=annotation_encoder,
    )

    expected_encoding = [  # TG**
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
    ]

    test_variant_snp = MockVariantRecord(
        contig="contig02",
        pos=1,
        ref="T",
        alts=("C",),
        info={
            "MQ": 3.1,
        },
    )
    observed_representation = reference_reader.variant_representation(
        test_variant_snp
    )
    assert observed_representation["chrom"] == "contig02"
    assert observed_representation["pos"] == 1
    assert observed_representation["ref"] == "T"
    assert observed_representation["alt"] == "C"
    assert observed_representation["type"] == 0.0
    npt.assert_array_almost_equal(
        observed_representation["best_practices"], [3.1, 0.0]
    )
    npt.assert_array_equal(
        observed_representation["reference"],
        expected_encoding,
    )

    test_variant_indel = MockVariantRecord(
        contig="contig02",
        pos=1,
        ref="T",
        alts=("C", "CT"),
        info={"MQ": 3.1, "DP": 1.4},
    )
    observed_representation = reference_reader.variant_representation(
        test_variant_indel
    )
    assert observed_representation["chrom"] == "contig02"
    assert observed_representation["pos"] == 1
    assert observed_representation["ref"] == "T"
    assert observed_representation["alt"] == "C, CT"
    assert observed_representation["type"] == -1.0
    npt.assert_array_almost_equal(
        observed_representation["best_practices"], [3.1, 1.4]
    )
    npt.assert_array_equal(
        observed_representation["reference"],
        expected_encoding,
    )

    observed_representation = reference_reader.variant_representation(
        test_variant_indel, alt_sep=" "
    )
    assert observed_representation["alt"] == "C CT"

    test_variant_indel = MockVariantRecord(
        contig="contig02",
        pos=1,
        ref="T",
        alts=("CG", "CT"),
        info={"MQ": 3.1, "DP": 1.4},
    )
    observed_representation = reference_reader.variant_representation(
        test_variant_indel, alt_sep=" "
    )
    assert observed_representation["type"] == 1.0


def test_reference_tensor_reader():
    annotation_encoder = AnnotationEncoder(annotation_list=["MQ", "DP"])

    reader = TensorReader(
        variant_file=None,
        reference_dict=seq_dict,
        read_encoder=lambda x: [[[7.25]]],
        annotation_encoder=annotation_encoder,
    )

    expected_encoding = [  # TG**
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
    ]

    test_variant_snp = MockVariantRecord(
        contig="contig02",
        pos=1,
        ref="T",
        alts=("C",),
        info={
            "MQ": 3.1,
        },
    )
    observed_representation = reader.variant_representation(
        test_variant_snp
    )
    assert observed_representation["chrom"] == "contig02"
    assert observed_representation["pos"] == 1
    assert observed_representation["ref"] == "T"
    assert observed_representation["alt"] == "C"
    assert observed_representation["type"] == 0.0
    npt.assert_array_almost_equal(
        observed_representation["best_practices"], [3.1, 0.0]
    )
    npt.assert_array_almost_equal(observed_representation["read_tensor"],
                                 [[[7.25]]])