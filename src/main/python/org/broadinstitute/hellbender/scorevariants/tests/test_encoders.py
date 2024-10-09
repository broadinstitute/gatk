from Bio import SeqIO
from scorevariants.encoders import (
    VariantLabelEncoder,
    VariantTypeEncoder,
    ReferenceEncoder,
    AnnotationEncoder,
    DNAVectorEncoder,
    ReadTensorEncoder,
    Interval,
    CIGAR_CODE,
)
from scorevariants.tests.utilities import (
    _make_test_fasta,
    MockVariantRecord,
)

import numpy as np
import numpy.testing as npt
import unittest

tc = unittest.TestCase()


def test_dna_vector_encoder():
    encoder = DNAVectorEncoder()
    obs = encoder("A")
    npt.assert_array_equal(obs, [1.0, 0, 0, 0])
    obs = encoder("N")
    npt.assert_array_equal(obs, [0.25, 0.25, 0.25, 0.25])


def test_variant_label_encoder():
    encoder = VariantLabelEncoder()
    assert encoder("SNP") == 2
    assert encoder("NOT_SNP") == 0
    assert encoder("INDEL") == 3
    assert encoder("NOT_INDEL") == 1
    assert encoder("NONE") == -1
    assert encoder.inverse(0) == "NOT_SNP"
    assert encoder.inverse(1) == "NOT_INDEL"
    assert encoder.inverse(2) == "SNP"
    assert encoder.inverse(3) == "INDEL"
    assert encoder.inverse(-1) == "NONE"


class MockAlignmentRead:
    def __init__(
        self,
        tags=None,
        reference_start=None,
        cigarstring=None,
        cigartuples=None,
        query_alignment_sequence=None,
        query_alignment_qualities=None,
    ):
        if cigarstring:
            self.cigarstring = cigarstring
        self.tags = tags
        self.reference_start = reference_start
        self.cigartuples = cigartuples
        self.query_alignment_sequence = query_alignment_sequence
        self.query_alignment_qualities = query_alignment_qualities

    def get_tag(self, key):
        return self.tags.get(key)


def test_variant_type_encoder():
    encoder = VariantTypeEncoder()
    assert encoder("SNP") == 0
    assert encoder("INDEL") == 1
    assert encoder("OTHER") == -1
    assert encoder.inverse(0) == "SNP"
    assert encoder.inverse(1) == "INDEL"
    assert encoder.inverse(-1) == "OTHER"


def test_read_tensor_encoder_init():
    encoder = ReadTensorEncoder(None, None)
    assert encoder.n_channels == 15
    channel_map = encoder.channel_map
    assert channel_map["read_A"] == 0
    assert channel_map["read_C"] == 1
    assert channel_map["read_G"] == 2
    assert channel_map["read_T"] == 3
    assert channel_map["read_*"] == 4
    assert channel_map["reference_A"] == 5
    assert channel_map["reference_C"] == 6
    assert channel_map["reference_G"] == 7
    assert channel_map["reference_T"] == 8
    assert channel_map["reference_*"] == 9
    assert channel_map["flag_bit_4"] == 10
    assert channel_map["flag_bit_4"] == 10
    assert channel_map["flag_bit_5"] == 11
    assert channel_map["flag_bit_6"] == 12
    assert channel_map["flag_bit_7"] == 13
    assert channel_map["mapping_quality"] == 14


def test_read_tensor_encoder_initialize_representation():
    encoder = ReadTensorEncoder(None, None)
    tensor = encoder.initialize_representation()
    assert np.all(tensor == 0.0)
    assert tensor.shape == (
        128,
        128,
        15,
    )


def test_read_tensor_encoder_read_filters():
    encoder = ReadTensorEncoder(None, None)
    interval = Interval(1000, 1128)
    mock_reads = [
        (
            MockAlignmentRead(reference_start=1111, tags={"RG": "some_group"}),
            True,
        ),
        (
            MockAlignmentRead(
                cigarstring=None,
                reference_start=1111,
                tags={"RG": "some_group"},
            ),
            True,
        ),
        (
            MockAlignmentRead(
                cigarstring="3M1I3M1D5M",
                reference_start=1111,
                tags={"RG": "some_group"},
            ),
            False,
        ),
        (
            MockAlignmentRead(
                cigarstring="3M1I3M1D5M",
                reference_start=1111,
                tags={"RG": "ArtificialHaplotype"},
            ),
            True,
        ),
        # left out for now while I resolve difference with gatk
        # (
        #     MockAlignmentRead(
        #         cigarstring='3M1I3M1D5M', reference_start=200, tags={'RG': 'some_group'}
        #      ),
        #     True
        # ),
    ]

    for read, exp in mock_reads:
        assert encoder.filter_read(read, interval) == exp


def test_read_tensor_encoder_get_insertions():
    encoder = ReadTensorEncoder(None, None)
    interval = Interval(1000, 1128)
    mock_reads = [
        MockAlignmentRead(
            cigarstring="3M1I3M1D5M",
            cigartuples=[
                (CIGAR_CODE["M"], 3),
                (CIGAR_CODE["I"], 1),
                (CIGAR_CODE["M"], 3),
                (CIGAR_CODE["D"], 1),
                (CIGAR_CODE["M"], 5),
            ],
            reference_start=1100,
        ),
        MockAlignmentRead(
            cigarstring="4M5I3M2I5M",
            cigartuples=[
                (CIGAR_CODE["M"], 4),
                (CIGAR_CODE["I"], 5),
                (CIGAR_CODE["M"], 3),
                (CIGAR_CODE["I"], 2),
                (CIGAR_CODE["M"], 5),
            ],
            reference_start=1004,
        ),
    ]
    expectations = [
        {103: 1},
        {
            8: 5,
            16: 2,
        },
    ]
    assert len(mock_reads) == len(expectations)

    for read, exp in zip(mock_reads, expectations):
        insertions = encoder.get_insertions(read, interval)
        tc.assertDictEqual(insertions, exp)


def test_read_tensor_encoder_get_insertions_back_to_back():
    encoder = ReadTensorEncoder(None, None)
    interval = Interval(1000, 1128)
    inserts = dict()
    mock_reads = [
        MockAlignmentRead(
            cigarstring="4M5I3M2I5M",
            cigartuples=[
                (CIGAR_CODE["M"], 4),
                (CIGAR_CODE["I"], 5),
                (CIGAR_CODE["M"], 3),
                (CIGAR_CODE["I"], 2),
                (CIGAR_CODE["M"], 5),
            ],
            reference_start=1004,
        ),
        MockAlignmentRead(
            cigarstring="4M4I3M2I5M",
            cigartuples=[
                (CIGAR_CODE["M"], 4),
                (CIGAR_CODE["I"], 4),
                (CIGAR_CODE["M"], 3),
                (CIGAR_CODE["I"], 2),
                (CIGAR_CODE["M"], 5),
            ],
            reference_start=1004,
        ),
    ]
    for read in mock_reads:
        encoder.get_insertions(read, interval, inserts)
    expectation = {
        8: 5,
        16: 2,
        15: 2,
    }
    tc.assertDictEqual(expectation, inserts)


def test_read_tensor_encoder_get_reference():
    reference = {"contig1": "ACGTGGGGTTTT"}
    encoder = ReadTensorEncoder(None, None)
    interval = Interval(1, 5)
    contig = "contig1"
    insertions = {
        -2: 2,
        2: 5,
    }
    exp = "**CG*****TG"
    reference_seq = encoder.get_reference(
        reference, contig, interval, insertions
    )

    assert exp == reference_seq
    
        
    interval = Interval(-2, 3)
    exp = 'ACG'
    reference_seq =encoder.get_reference(
           reference, contig, interval, {} 
    )
    assert exp == reference_seq


def test_read_tensor_encoder_get_sequence_specifc_indel_dict():
    encoder = ReadTensorEncoder(None, None)
    # interval = Interval(1000, 1128)
    offset = -4
    all_read_insertions = {
        8: 5,
        16: 2,
    }

    # Total match to all read insertions
    read = MockAlignmentRead(
        cigarstring="4M5I3M2I5M",
        cigartuples=[
            (CIGAR_CODE["M"], 4),
            (CIGAR_CODE["I"], 5),
            (CIGAR_CODE["M"], 3),
            (CIGAR_CODE["I"], 2),
            (CIGAR_CODE["M"], 5),
        ],
        reference_start=1004,
    )
    observed_indels = encoder._get_sequence_specific_indel_dict(
        read, offset, all_read_insertions
    )
    # the read has all insertions from the reads,
    # so the read does not need *'s inserted
    assert sum(abs(v) for v in observed_indels.values()) == 0

    # no match to all read insertions
    read = MockAlignmentRead(
        cigarstring="19M",
        cigartuples=[
            (CIGAR_CODE["M"], 19),
        ],
        reference_start=1004,
    )
    observed_indels = encoder._get_sequence_specific_indel_dict(
        read, offset, all_read_insertions
    )
    # the read matches the reference, so it needs all of the *'s
    # to be inserted
    tc.assertDictEqual(observed_indels, all_read_insertions)

    # partial match of insertions
    read = MockAlignmentRead(
        cigarstring="4M4I3M2I5M",
        cigartuples=[
            (CIGAR_CODE["M"], 4),
            (CIGAR_CODE["I"], 4),
            (CIGAR_CODE["M"], 3),
            (CIGAR_CODE["I"], 2),
            (CIGAR_CODE["M"], 5),
        ],
        reference_start=1004,
    )
    exp_insertions = {
        8: 1,
        15: 0,
        16: 2,
    }
    all_read_insertions_2 = {
        8: 5,
        15: 2,
        16: 2,
    }
    observed_indels = encoder._get_sequence_specific_indel_dict(
        read, offset, all_read_insertions_2
    )
    tc.assertDictEqual(exp_insertions, observed_indels)

    # adjusts for deletions
    read = MockAlignmentRead(
        cigarstring="4M5I1D2M2I5M",
        cigartuples=[
            (CIGAR_CODE["M"], 4),
            (CIGAR_CODE["I"], 5),
            (CIGAR_CODE["D"], 1),
            (CIGAR_CODE["M"], 2),
            (CIGAR_CODE["I"], 1),
            (CIGAR_CODE["M"], 5),
        ],
        reference_start=1004,
    )
    exp_insertions = {
        8: 0,
        13: 1,
        16: 1,
    }
    observed_indels = encoder._get_sequence_specific_indel_dict(
        read, offset, all_read_insertions
    )
    tc.assertDictEqual(exp_insertions, observed_indels)


def test_read_tensor_encoder_get_sequence_and_quality():
    encoder = ReadTensorEncoder(None, None, window=10)
    interval = Interval(1000, 1010)
    read = MockAlignmentRead(
        cigartuples=[
            (CIGAR_CODE["M"], 2),
            (CIGAR_CODE["I"], 2),
            (CIGAR_CODE["M"], 2),
        ],
        reference_start=1004,
        query_alignment_sequence="ACGTACG",
        query_alignment_qualities=np.array(
            [
                7,
                8,
                9,
                10,
                11,
                12,
                13,
            ]
        ),
    )
    insertions = {
        6: 3,
    }
    seq1, qual1 = encoder._get_sequence_and_quality(read, insertions, interval)

    assert seq1[:10] == "~~~~AC*GTA"
    assert qual1[:10] == [0, 0, 0, 0, 7, 8, 0, 9, 10, 11]

    read = MockAlignmentRead(
        cigartuples=[
            (CIGAR_CODE["M"], 3),
        ],
        reference_start=1003,
        query_alignment_sequence="TAC",
        query_alignment_qualities=np.array(
            [
                7,
                8,
                9,
            ]
        ),
    )
    insertions = {
        5: 3,
    }
    seq2, qual2 = encoder._get_sequence_and_quality(read, insertions, interval)

    assert seq2 == "~~~TA***C"
    assert qual2 == [0, 0, 0, 7, 8, 0, 0, 0, 9]


def test_read_tensor_encoder_insert_reference():
    encoder = ReadTensorEncoder(None, None, read_limit=2, window=4)
    encoder.channel_map = {
        f"reference_{base}": i for i, base in enumerate("ACGT*")
    }

    tensor = np.zeros((2, 4, 5))
    reference = "AC*T"
    expectation = 2 * [
        [
            [1.0, 0, 0, 0, 0],
            [0.0, 1, 0, 0, 0],
            [0.0, 0, 0, 0, 1],
            [0.0, 0, 0, 1, 0],
        ]
    ]
    encoder.insert_reference(reference, tensor)
    npt.assert_array_equal(expectation, tensor)

    tensor = np.zeros((2, 4, 5))
    reference = "**KC*T"
    expectation = 2 * [
        [
            [0.0, 0, 0, 0, 1],
            [0.0, 0, 0, 0, 1],
            [0.0, 0, 0.5, 0.5, 0],
            [0.0, 1, 0, 0, 0],
        ]
    ]
    encoder.insert_reference(reference, tensor)
    npt.assert_array_equal(expectation, tensor)
    
    tensor = np.zeros((2, 4, 5))
    reference = "**KC*T"
    expectation = 2 * [
        [
            [0.0, 0, 0, 0, 0],
            [0.0, 0, 0, 0, 1],
            [0.0, 0, 0, 0, 1],
            [0.0, 0, 0.5, 0.5, 0],
        ]
    ]
    encoder.insert_reference(reference, tensor, start=1)
    npt.assert_array_equal(expectation, tensor)


def test_annotation_encoder():
    annotation_encoder = AnnotationEncoder(annotation_list=["MQ", "DP"])
    mock_variants = [
        MockVariantRecord(
            1,
            "c1",
            {
                "DP": 1.0,
                "MQ": 2,
            },
        ),
        MockVariantRecord(
            2,
            "c2",
            {
                "DP": 3.33,
            },
        ),
    ]
    encodings = [annotation_encoder(variant) for variant in mock_variants]
    expected_encodings = [
        [2.0, 1.0],
        [0.0, 3.33],
    ]
    npt.assert_array_equal(expected_encodings, encodings)


def test_reference_encoder():
    ref_encoder = ReferenceEncoder(window=6)
    reference = SeqIO.to_dict(
        SeqIO.parse(
            _make_test_fasta(),
            format="fasta",
        )
    )
    mock_variant = MockVariantRecord(pos=4, contig="SEQUENCE_1")
    expected_encoding = [  # ACGTTC
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0, 0.0],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)


def test_reference_encoder_with_offset_and_off_right_side():
    ref_encoder = ReferenceEncoder(window=6, offset=2)
    reference = SeqIO.to_dict(
        SeqIO.parse(
            _make_test_fasta(),
            format="fasta",
        )
    )
    mock_variant = MockVariantRecord(pos=4, contig="SEQUENCE_1")
    expected_encoding = [  # GTTCG*
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)


def test_reference_encoder_off_left_side():
    ref_encoder = ReferenceEncoder(window=4, handle_start="pad")
    reference = SeqIO.to_dict(
        SeqIO.parse(
            _make_test_fasta(),
            format="fasta",
        )
    )
    mock_variant = MockVariantRecord(pos=1, contig="SEQUENCE_2")
    expected_encoding = [  # **TG
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)

    ref_encoder = ReferenceEncoder(window=4)  # default to handle_start='left'
    reference = SeqIO.to_dict(
        SeqIO.parse(
            _make_test_fasta(),
            format="fasta",
        )
    )
    mock_variant = MockVariantRecord(pos=1, contig="SEQUENCE_2")
    expected_encoding = [  # TGC*
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)


def test_reference_encoder_ambiguous_encoding():
    ref_encoder = ReferenceEncoder(window=4)  # default to handle_start='left'
    reference = SeqIO.to_dict(
        SeqIO.parse(
            _make_test_fasta(),
            format="fasta",
        )
    )
    mock_variant = MockVariantRecord(pos=1, contig="SEQUENCE_3")
    expected_encoding = [  # NGC*
        [0.25, 0.25, 0.25, 0.25],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)

    mock_variant = MockVariantRecord(pos=2, contig="SEQUENCE_3")
    expected_encoding = [  # NGCR
        [0.25, 0.25, 0.25, 0.25],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.25, 0.25, 0.25, 0.25],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)

    ref_encoder = ReferenceEncoder(
        window=4, base_encoder=DNAVectorEncoder("code")
    )  # default to handle_start='left'
    mock_variant = MockVariantRecord(pos=2, contig="SEQUENCE_3")
    expected_encoding = [  # NGCR
        [0.25, 0.25, 0.25, 0.25],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.5, 0.0, 0.0, 0.5],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)


def test_reference_encoder_off_right_side():
    ref_encoder = ReferenceEncoder(window=4)
    reference = SeqIO.to_dict(
        SeqIO.parse(
            _make_test_fasta(),
            format="fasta",
        )
    )
    mock_variant = MockVariantRecord(pos=7, contig="SEQUENCE_1")
    expected_encoding = [  # TCG*
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0],
    ]
    encoding = ref_encoder(reference, mock_variant)
    npt.assert_array_equal(expected_encoding, encoding)
