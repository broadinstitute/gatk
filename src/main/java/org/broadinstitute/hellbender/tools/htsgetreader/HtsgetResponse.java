package org.broadinstitute.hellbender.tools.htsgetreader;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Base64;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.fasterxml.jackson.annotation.JsonProperty;
import com.fasterxml.jackson.annotation.JsonRootName;
import com.fasterxml.jackson.databind.annotation.JsonDeserialize;

import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.HttpUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

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
                    try (final CloseableHttpResponse resp = HttpUtils.getClient().execute(get)) {
                        final Path outputFile = IOUtils.createTempPath("htsget-temp", "");
                        try (final OutputStream ostream = Files.newOutputStream(outputFile);
                             final InputStream istream = resp.getEntity().getContent()) {
                            org.apache.commons.io.IOUtils.copy(istream, ostream);
                        } catch (final IOException e) {
                            throw new UserException("Could not write to temp file", e);
                        }
                        return Files.newInputStream(outputFile);
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
}
