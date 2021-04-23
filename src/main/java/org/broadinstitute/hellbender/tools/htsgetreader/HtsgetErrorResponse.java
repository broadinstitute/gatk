package org.broadinstitute.hellbender.tools.htsgetreader;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonRootName;

/**
 * Class allowing deserialization from json htsget error response
 */
@JsonRootName(value = "htsget")
public class HtsgetErrorResponse {
    @JsonProperty("error")
    private String error;

    @JsonProperty("message")
    private String message;

    public String getError() {
        return error;
    }

    public String getMessage() {
        return message;
    }
}
