package org.broadinstitute.hellbender.tools.htsgetreader;

import java.io.*;
import java.net.*;
import java.util.ArrayList;
import java.util.Base64;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Stream;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonRootName;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;

import org.apache.commons.io.input.AutoCloseInputStream;
import org.apache.http.HttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.HttpUtils;

/**
 * Class allowing deserialization from json htsget response
 */
@JsonRootName(value = "htsget")
public class HtsgetResponse {
    public static class Block {
        @JsonProperty("url")
        private URI uri;

        @JsonProperty("headers")
        @JsonDeserialize(as = HashMap.class, keyAs = String.class, contentAs = String.class)
        private Map<String, String> headers;

        @JsonProperty("class")
        private HtsgetClass dataClass;

        public URI getUri() {
            return this.uri;
        }

        public Map<String, String> getHeaders() {
            return Collections.unmodifiableMap(this.headers);
        }

        public HtsgetClass getDataClass() {
            return this.dataClass;
        }

        /**
         * Returns InputStream containing data of a single block from the htsget response
         *
         * Large blocks (those behind an http(s) URI) are first saved
         * to a temp file that is deleted upon program exit
         */
        public InputStream getData() {
            switch (this.getUri().getScheme()) {
                case "http":
                case "https":
                    final HttpGet get = new HttpGet(this.getUri());
                    this.getHeaders().forEach(get::addHeader);
                    try {
                        final HttpResponse resp = HttpUtils.getClient().execute(get);
                        return new AutoCloseInputStream(resp.getEntity().getContent());
                    } catch (final IOException e) {
                        throw new UserException("Could not retrieve data from block", e);
                    }
                case "data":
                    final String dataUri = this.getUri().toString();
                    if (!dataUri.matches("^data:.*;base64,.*")) {
                        throw new UserException("data URI must be base64 encoded: " + dataUri);
                    }
                    return new ByteArrayInputStream(
                            Base64.getDecoder().decode(dataUri.replaceFirst("^data:.*;base64,", "")));
                default:
                    throw new UserException("Unrecognized URI scheme in data block: " + this.getUri().getScheme());
            }
        }
    }

    @JsonProperty("format")
    private HtsgetFormat format;

    @JsonProperty("urls")
    @JsonDeserialize(as = ArrayList.class, contentAs = Block.class)
    private List<Block> blocks;

    @JsonProperty("md5")
    private String md5;

    public HtsgetFormat getFormat() {
        return this.format;
    }

    public List<Block> getBlocks() {
        return Collections.unmodifiableList(this.blocks);
    }

    public String getMd5() {
        return this.md5;
    }

    public Stream<InputStream> streamData() {
        return this.getBlocks().stream().map(HtsgetResponse.Block::getData);
    }
}
