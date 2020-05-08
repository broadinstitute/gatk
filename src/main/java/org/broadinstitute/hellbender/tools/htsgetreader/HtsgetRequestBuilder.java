package org.broadinstitute.hellbender.tools.htsgetreader;

import htsjdk.samtools.util.FileExtensions;

import java.net.URI;
import java.util.Collection;
import java.util.Collections;
import java.util.EnumSet;
import java.util.HashSet;
import java.util.Set;
import java.util.stream.Collectors;

import javax.ws.rs.core.UriBuilder;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;

/**
 * Builder for an htsget request that allows converting the request
 * to a URI after validating that it is properly formed
 */
public class HtsgetRequestBuilder {
    final private URI endpoint;
    final private String id;

    // Query parameters
    private HtsgetFormat format;
    private HtsgetClass dataClass;
    private SimpleInterval interval;
    private final EnumSet<HtsgetRequestField> fields;
    private final Set<String> tags;
    private final Set<String> notags;

    public HtsgetRequestBuilder(final URI endpoint, final String id) {
        this.endpoint = endpoint;
        this.id = id;
        this.fields = EnumSet.noneOf(HtsgetRequestField.class);
        this.tags = new HashSet<>();
        this.notags = new HashSet<>();
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

    public HtsgetRequestBuilder withFormat(final HtsgetFormat format) {
        this.format = format;
        return this;
    }

    public HtsgetRequestBuilder withDataClass(final HtsgetClass dataClass) {
        this.dataClass = dataClass;
        return this;
    }

    public HtsgetRequestBuilder withInterval(final SimpleInterval interval) {
        this.interval = interval;
        return this;
    }

    public HtsgetRequestBuilder withField(final HtsgetRequestField field) {
        this.fields.add(field);
        return this;
    }

    public HtsgetRequestBuilder withFields(final Collection<HtsgetRequestField> fields) {
        this.fields.addAll(fields);
        return this;
    }

    public HtsgetRequestBuilder withTag(final String tag) {
        this.tags.add(tag);
        return this;
    }

    public HtsgetRequestBuilder withTags(final Collection<String> tags) {
        this.tags.addAll(tags);
        return this;
    }

    public HtsgetRequestBuilder withNotag(final String notag) {
        this.notags.add(notag);
        return this;
    }

    public HtsgetRequestBuilder withNotags(final Collection<String> notags) {
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
            builder.queryParam("start", this.interval.getGA4GHStart());
            builder.queryParam("end", this.interval.getGA4GHEnd());
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
}
