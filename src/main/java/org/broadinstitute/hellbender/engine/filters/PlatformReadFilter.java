package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import org.broadinstitute.hellbender.cmdline.Argument;

import java.util.Set;

/**
 * Keep only reads that match th PL attribute.
 * Matching is done by case-insensitive substring matching
 * (checking if the read's platform tag contains the given string).
 */
public final class PlatformReadFilter implements ReadFilter {
    private static final long serialVersionUID = 1L;
    @Argument(fullName = "PLFilterName", shortName = "PLFilterName", doc="Keep reads with RG:PL attribute containing this string", optional=true)
    public Set<String> PLFilterNames;

    @Override
    public boolean test(final SAMRecord read) {
        final SAMReadGroupRecord readGroup = read.getReadGroup();
        if (readGroup == null) {
            return false;
        }
        final String readPlatformAttr = readGroup.getPlatform();
        if (readPlatformAttr == null) {
            return false;
        }
        final String platformUppercase = readPlatformAttr.toUpperCase();
        for (final String name : PLFilterNames) {
            if (platformUppercase.contains(name.toUpperCase())){
                return true;
            }
        }
        return false;
    }
}
