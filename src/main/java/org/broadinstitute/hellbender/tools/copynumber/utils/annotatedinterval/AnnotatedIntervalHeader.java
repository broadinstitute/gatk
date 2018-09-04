package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedinterval;

import htsjdk.samtools.SAMFileHeader;

import java.util.List;

public class AnnotatedIntervalHeader {
    private final String contigColumnName;
    private final String startColumnName;
    private final String endColumnName;
    private final List<String> annotations;
    private final SAMFileHeader samFileHeader;

    /**
     * @param samFileHeader SAM file header as a multiline string.  {@code null} is allowed, if not available.
     * @param annotations annotation names that do not include the locatable column names.  Never {@code null}.
     * @param contigColumnName how contig should be rendered.  Never {@code null}.
     * @param startColumnName how start position should be rendered.  Never {@code null}.
     * @param endColumnName how end position should be rendered.  Never {@code null}.
     */
    public AnnotatedIntervalHeader(final String contigColumnName, final String startColumnName, final String endColumnName, final List<String> annotations, final SAMFileHeader samFileHeader) {
        this.contigColumnName = contigColumnName;
        this.startColumnName = startColumnName;
        this.endColumnName = endColumnName;
        this.annotations = annotations;
        this.samFileHeader = samFileHeader;
    }


    public String getContigColumnName() {
        return contigColumnName;
    }

    public String getStartColumnName() {
        return startColumnName;
    }

    public String getEndColumnName() {
        return endColumnName;
    }

    public List<String> getAnnotations() {
        return annotations;
    }

    /** Can be {@code null} */
    public SAMFileHeader getSamFileHeader() {
        return samFileHeader;
    }
}
