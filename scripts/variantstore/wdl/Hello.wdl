version 1.0

import "GvsUtils.wdl" as Utils

workflow hello {
    call Utils.GetToolVersions
}