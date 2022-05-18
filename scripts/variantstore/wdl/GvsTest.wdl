version 1.0

import "https://raw.githubusercontent.com/broadinstitute/gatk/41943ca6634f888b765d4455e95df1ba916ee5d2/scripts/variantstore/wdl/GvsUtils.wdl" as GvsUtil

workflow GvsTest {
    call GvsUtil.TerminateWorkflow {
        input: message = "Holy moly 'here' are a bunch of \"problematic\" quotes"
    }
}
