package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

import com.google.common.annotations.VisibleForTesting;

@VisibleForTesting
public enum TypeInferredFromSimpleChimera {
    SIMPLE_DEL,             // simple deletion
    DEL_DUP_CONTRACTION,    // deletion from duplication contraction
    SIMPLE_INS,             // simple insertion
    RPL,                    // replacement
    SMALL_DUP_EXPANSION,    // small duplication tandem expansion
    SMALL_DUP_CPX,          // small duplication, either expansion or contraction, with complex structure

    INTRA_CHR_STRAND_SWITCH_55,// intra-chromosome strand-switch novel adjacency, alignments left-flanking the novel adjacency
    INTRA_CHR_STRAND_SWITCH_33,// intra-chromosome strand-switch novel adjacency, alignments right-flanking the novel adjacency

    INTRA_CHR_REF_ORDER_SWAP,// intra-chromosome reference-order swap, but NO strand-switch, novel adjacency

    INTER_CHR_STRAND_SWITCH_55,// pair WY in Fig.1 in Section 5.4 of VCF spec ver.4.2
    INTER_CHR_STRAND_SWITCH_33,// pair XZ in Fig.1 in Section 5.4 of VCF spec ver.4.2
    INTER_CHR_NO_SS_WITH_LEFT_MATE_FIRST_IN_PARTNER, // the green pair in Fig. 7 in Section 5.4 of VCF spec ver.4.2
    INTER_CHR_NO_SS_WITH_LEFT_MATE_SECOND_IN_PARTNER; // the red pair in Fig. 7 in Section 5.4 of VCF spec ver.4.2
}
