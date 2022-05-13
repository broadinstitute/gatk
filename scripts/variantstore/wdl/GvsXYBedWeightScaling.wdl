version 1.0

import "GvsUtils.wdl" as GvsUtils

workflow GvsXYBedWeightScaling {
    input {
        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
        Float scale_factor = 4
    }

    call GvsUtils.ScaleXYBedValues {
        input:
            interval_weights_bed = interval_weights_bed,
            scale_factor = scale_factor
    }
}
