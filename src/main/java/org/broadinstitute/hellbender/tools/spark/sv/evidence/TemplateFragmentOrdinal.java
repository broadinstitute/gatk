package org.broadinstitute.hellbender.tools.spark.sv.evidence;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.GATKRead;

/**
 * Indicates the ordinal of a fragment in a paired sequenced template.
 */
public enum TemplateFragmentOrdinal {

    /**
     * Ordinal for fragments in unpaired sequenced templates.
     */
    UNPAIRED(""),

    /**
     * For fragment in paired sequenced templates with unknown ordinal.
     */
    PAIRED_UNKNOWN("/?"),


    /**
     * For the first and only first fragment in a paired sequenced template.
     */
    PAIRED_FIRST("/1"),

    /**
     * For the second/last and only second/last fragment in a paired sequenced template.
     */
    PAIRED_SECOND("/2"),


    /**
     * For inner fragments in a paired sequenced template.
     */
    PAIRED_INTERIOR("/0");

    TemplateFragmentOrdinal(final String nameSuffix) {
        this.nameSuffix = nameSuffix;
    }

    @Override
    public String toString() {
        return nameSuffix;
    }

    public String nameSuffix() {
        return nameSuffix;
    }

    private final String nameSuffix;

    public static TemplateFragmentOrdinal forRead(final GATKRead read) {
        Utils.nonNull(read);
        if (read.isPaired()) {
            if (read.isFirstOfPair()) {
                return read.isSecondOfPair() ? PAIRED_INTERIOR : PAIRED_FIRST;
            } else {
                return read.isSecondOfPair() ? PAIRED_SECOND : PAIRED_UNKNOWN;
            }
        } else {
            return UNPAIRED;
        }
    }
}
