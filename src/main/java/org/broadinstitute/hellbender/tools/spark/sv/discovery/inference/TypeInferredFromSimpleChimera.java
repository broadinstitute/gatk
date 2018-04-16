package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

enum TypeInferredFromSimpleChimera {
    SIMPLE_DEL,             // simple deletion
    DEL_DUP_CONTRACTION,    // deletion from duplication contraction
    SIMPLE_INS,             // simple insertion
    RPL,                    // replacement
    SMALL_DUP_EXPANSION,    // small duplication tandem expansion
    SMALL_DUP_CPX,          // small duplication, either expansion or contraction, with complex structure
    IntraChrStrandSwitch,   // intra-chromosome strand-switch novel adjacency
    IntraChrRefOrderSwap,   // intra-chromosome reference-order swap, but NO strand-switch, novel adjacency
    InterChromosome;        // inter-chromosome novel adjacency (regardless of with or without strand-switch)
}
