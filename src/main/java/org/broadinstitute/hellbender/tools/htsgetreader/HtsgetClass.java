package org.broadinstitute.hellbender.tools.htsgetreader;

import com.fasterxml.jackson.annotation.JsonProperty;

/**
 * Classes of data that can be requested in an htsget request as defined by the spec
 */
public enum HtsgetClass {
    @JsonProperty("body")
    body,

    @JsonProperty("header")
    header;
}
