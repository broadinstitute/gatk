package org.broadinstitute.hellbender.tools.htsgetreader;

import com.fasterxml.jackson.databind.DeserializationFeature;
import com.fasterxml.jackson.databind.MapperFeature;
import com.fasterxml.jackson.databind.ObjectMapper;
import htsjdk.samtools.util.FileExtensions;

import java.io.IOException;
import java.net.URI;
import java.net.URISyntaxException;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import javax.ws.rs.core.UriBuilder;

import org.apache.commons.io.Charsets;
import org.apache.http.Header;
import org.apache.http.HttpEntity;
import org.apache.http.client.methods.CloseableHttpResponse;
import org.apache.http.client.methods.HttpGet;
import org.apache.http.util.EntityUtils;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.HttpUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * Builder for an htsget request that allows converting the request
 * to a URI after validating that it is properly formed
 */
public class HtsgetRequest {
    final private URI endpoint;
    final private String id;

    // Query parameters
    private HtsgetFormat format;
    private HtsgetClass dataClass;
    private SimpleInterval interval;
    private final EnumSet<HtsgetRequestField> fields;
    private final Set<String> tags;
    private final Set<String> notags;

    public HtsgetRequest(final URI endpoint, final String id) {
        this.endpoint = endpoint;
        this.id = id;
        this.fields = EnumSet.noneOf(HtsgetRequestField.class);
        this.tags = new HashSet<>();
        this.notags = new HashSet<>();
    }

    public HtsgetRequest(final GATKPath source) {
        try {
            final URI sourceURI = source.getURI();
            this.endpoint = new URI("//" + sourceURI.getHost());
            this.id = sourceURI.getPath();
            this.fields = EnumSet.noneOf(HtsgetRequestField.class);
            this.tags = new HashSet<>();
            this.notags = new HashSet<>();
        } catch (final URISyntaxException e) {
            throw new UserException(source.toString(), e);
        }
    }

    public URI getEndpoint() {
        return this.endpoint;
    }

    public String getID() {
        return this.id;
    }

    public HtsgetFormat getFormat() {
        return this.format;
    }

    public HtsgetClass getDataClass() {
        return this.dataClass;
    }

    public SimpleInterval getInterval() {
        return this.interval;
    }

    public Set<HtsgetRequestField> getFields() {
        return Collections.unmodifiableSet(this.fields);
    }

    public Set<String> getTags() {
        return Collections.unmodifiableSet(this.tags);
    }

    public Set<String> getNoTags() {
        return Collections.unmodifiableSet(this.notags);
    }

    public void setFormat(final HtsgetFormat format) {
        this.format = format;
    }

    public void setDataClass(final HtsgetClass dataClass) {
        this.dataClass = dataClass;
    }

    public void setInterval(final SimpleInterval interval) {
        this.interval = interval;
    }

    public void addField(final HtsgetRequestField field) {
        this.fields.add(field);
    }

    public void addFields(final Collection<HtsgetRequestField> fields) {
        this.fields.addAll(fields);
    }

    public void addTag(final String tag) {
        this.tags.add(tag);
    }

    public void addTags(final Collection<String> tags) {
        this.tags.addAll(tags);
    }

    public void addNotag(final String notag) {
        this.notags.add(notag);
    }

    public void addNotags(final Collection<String> notags) {
        this.notags.addAll(notags);
    }

    public HtsgetRequest withFormat(final HtsgetFormat format) {
        this.format = format;
        return this;
    }

    public HtsgetRequest withDataClass(final HtsgetClass dataClass) {
        this.dataClass = dataClass;
        return this;
    }

    public HtsgetRequest withInterval(final SimpleInterval interval) {
        this.interval = interval;
        return this;
    }

    public HtsgetRequest withField(final HtsgetRequestField field) {
        this.fields.add(field);
        return this;
    }

    public HtsgetRequest withFields(final Collection<HtsgetRequestField> fields) {
        this.fields.addAll(fields);
        return this;
    }

    public HtsgetRequest withTag(final String tag) {
        this.tags.add(tag);
        return this;
    }

    public HtsgetRequest withTags(final Collection<String> tags) {
        this.tags.addAll(tags);
        return this;
    }

    public HtsgetRequest withNotag(final String notag) {
        this.notags.add(notag);
        return this;
    }

    public HtsgetRequest withNotags(final Collection<String> notags) {
        this.notags.addAll(notags);
        return this;
    }

    /**
     * Validates that the user query obeys htsget spec
     */
    private void validateRequest() {
        if (this.dataClass != null && this.dataClass == HtsgetClass.header && (
            this.interval != null ||
            ! this.fields.isEmpty() ||
            ! this.tags.isEmpty() ||
            ! this.notags.isEmpty())) {
                throw new UserException("Invalid request: no query parameters except `format` may be specified when class=header");
        }

        if (this.format != null) {
            if ((this.id.endsWith(FileExtensions.BAM) || this.id.endsWith(FileExtensions.SAM) && (
                this.format != HtsgetFormat.BAM && this.format != HtsgetFormat.CRAM))
                ||
                FileExtensions.VCF_LIST.stream().anyMatch(this.id::endsWith) && (
                this.format != HtsgetFormat.VCF && this.format != HtsgetFormat.BCF)) {
                throw new UserException("Specified format: " + this.format + " is incompatible with id's file extension");
            }
        }

        final String intersections = this.tags.stream()
            .filter(getNoTags()::contains)
            .collect(Collectors.joining(", "));
        if (! intersections.isEmpty()) {
            throw new UserException("Invalid request: tags and notags overlap in the following fields: " + intersections);
        }
    }

    /**
     * Convert request to a URI which can be used to make http request for data blocks
     */
    public URI toURI() {
        this.validateRequest();
        final UriBuilder builder = UriBuilder.fromUri(this.endpoint)
            .scheme("https")
            .path(this.id);

        if (this.format != null) {
            builder.queryParam("format", this.format.toString());
        }
        if (this.dataClass != null) {
            builder.queryParam("class", this.dataClass);
        }
        if (this.interval != null) {
            builder.queryParam("referenceName", this.interval.getContig());
            // do not insert start and end for unmapped reads
            if (!this.interval.getContig().equals("*")) {
                builder.queryParam("start", this.interval.getGA4GHStart());
                builder.queryParam("end", this.interval.getGA4GHEnd());
            }
        }
        if (!this.fields.isEmpty()) {
            builder.queryParam(
                "fields",
                this.fields.stream().map(HtsgetRequestField::toString).collect(Collectors.joining(",")));
        }
        if (!this.tags.isEmpty()) {
            builder.queryParam("tags", String.join(",", this.tags));
        }
        if (!this.notags.isEmpty()) {
            builder.queryParam("notags", String.join(",", this.notags));
        }
        return builder.build();
    }

    private ObjectMapper getObjectMapper() {
        final ObjectMapper mapper = new ObjectMapper();
        mapper.enable(DeserializationFeature.UNWRAP_ROOT_VALUE);
        mapper.configure(MapperFeature.ACCEPT_CASE_INSENSITIVE_PROPERTIES, true);
        return mapper;
    }

    /**
     * Attempt to make htsget request and return response if there are no errors
     * @return the response from the htsget server if request is successful as an HtsgetResponse object
     */
    public HtsgetResponse getResponse() {
        final URI reqURI = this.toURI();

        final HttpGet getReq = new HttpGet(reqURI);
        try (final CloseableHttpResponse resp = HttpUtils.getClient().execute(getReq)) {
            // get content of response
            final HttpEntity entity = resp.getEntity();
            final Header encodingHeader = entity.getContentEncoding();
            final Charset encoding = encodingHeader == null
                ? StandardCharsets.UTF_8
                : Charsets.toCharset(encodingHeader.getValue());
            final String jsonBody = EntityUtils.toString(entity, encoding);

            final ObjectMapper mapper = this.getObjectMapper();

            if (resp.getStatusLine() == null) {
                throw new UserException("htsget server response did not contain status line");
            }
            final int statusCode = resp.getStatusLine().getStatusCode();
            if (400 <= statusCode && statusCode < 500) {
                final HtsgetErrorResponse err = mapper.readValue(jsonBody, HtsgetErrorResponse.class);
                throw new UserException(
                    "Invalid request, received error code: " + statusCode +
                        ", error type: " + err.getError() +
                        ", message: " + err.getMessage());
            } else if (statusCode == 200) {
                return mapper.readValue(jsonBody, HtsgetResponse.class);
            } else {
                throw new UserException("Unrecognized status code: " + statusCode);
            }
        } catch (final IOException e) {
            throw new UserException("IOException while attempting htsget download", e);
        }
    }
}
