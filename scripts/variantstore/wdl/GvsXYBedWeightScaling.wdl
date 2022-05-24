version 1.0

import "GvsUtils.wdl" as GvsUtils

workflow GvsXYBedWeightScaling {
    input {
        File interval_weights_bed = "gs://broad-public-datasets/gvs/weights/gvs_vet_weights_1kb.bed"
        Float x_bed_weight_scaling = 4.0
        Float y_bed_weight_scaling = 4.0
    }

    call GvsUtils.ScaleXYBedValues {
        input:
            interval_weights_bed = interval_weights_bed,
            x_bed_weight_scaling = x_bed_weight_scaling,
            y_bed_weight_scaling = y_bed_weight_scaling
    }
}
