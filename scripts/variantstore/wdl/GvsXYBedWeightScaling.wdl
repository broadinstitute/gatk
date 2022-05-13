version 1.0

import "GvsUtils.wdl" as GvsUtils

workflow GvsXYBedWeightScaling {
    call GvsUtils.ScaleXYBedValues {
        input:
            interval_weights_bed = "",
            scale_factor = 5
    }
}
