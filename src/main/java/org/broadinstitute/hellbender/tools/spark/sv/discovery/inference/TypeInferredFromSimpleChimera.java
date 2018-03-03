package org.broadinstitute.hellbender.tools.spark.sv.discovery.inference;

enum TypeInferredFromSimpleChimera {
    SIMPLE_DEL, DEL_DUP_CONTRACTION, SIMPLE_INS, RPL, SMALL_DUP_EXPANSION, SMALL_DUP_CPX,
    IntraChrStrandSwitch, IntraChrRefOrderSwap, InterChromosome;
}
