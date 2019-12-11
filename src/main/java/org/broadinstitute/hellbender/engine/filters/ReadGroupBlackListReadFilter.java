package org.broadinstitute.hellbender.engine.filters;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMTag;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.ReadFilterArgumentDefinitions;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.*;

/**
 * Keep records that don't match the specified filter string(s). Filter strings consist of a two character
 * read group tag name such as "RG", "PU", etc., as defined by {@link SAMReadGroupRecord}, followed by a
 * ":", and then the specific value to use for filtering.
 *
 * <p>For example, this filter value uses the platform unit (PU) tag:
 *   <code>PU:1000G-mpimg-080821-1_1</code>
 * to filter out reads with the read group platform unit value <code>1000G-mpimg-080821-1_1</code></p>
 */
@DocumentedFeature(groupName= HelpConstants.DOC_CAT_READFILTERS, groupSummary=HelpConstants.DOC_CAT_READFILTERS_SUMMARY)
public final class ReadGroupBlackListReadFilter extends ReadFilter implements Serializable {
    private static final long serialVersionUID = 1L;
    public static final String FILTER_ENTRY_SEPARATOR = ":";

    @Argument(fullName=ReadFilterArgumentDefinitions.READ_GROUP_BLACK_LIST_LONG_NAME,
            doc="A read group filter expression in the form \"attribute:value\", where \"attribute\" is "
                    + "a two character read group attribute such as \"RG\" or \"PU\".",
            optional=false)
    public List<String> blackList = new ArrayList<>();

    //most of the collection Entry classes are not serializable so just use a Map
    private final Map<String, Collection<String>> blacklistEntries = new HashMap<>();

    // Command line parser requires a no-arg constructor
    public ReadGroupBlackListReadFilter() {};

    @Override
    public void setHeader(final SAMFileHeader samFileHeader) {
        super.setHeader(samFileHeader);
        parseReadGroupFilters();
    }

    /**
     * Creates a filter using the list of blacklisted read groups.
     */
    public ReadGroupBlackListReadFilter(final List<String> blackLists, final SAMFileHeader header) {
        super.setHeader(header);
        this.blackList.addAll(blackLists);
        parseReadGroupFilters();
    }

    private void parseReadGroupFilters() {
        final Map<String, Collection<String>> filters = new TreeMap<>();
        for (String blackList : this.blackList) {
            addFilters(filters, blackList);
        }
        //merge all the new entries in to the blacklist
        filters.forEach((k, v) -> blacklistEntries.merge(k, v, (v1, v2) -> {
            v1.addAll(v2);
            return v1;
        }));
    }

    private void addFilters(final Map<String, Collection<String>> filters, final String filter) {
        final String[] split = filter.split(FILTER_ENTRY_SEPARATOR, 2);
        checkValidFilterEntry(filter, split);

        //Note: if we're here, we know that split has exactly 2 elements.
        filters.computeIfAbsent(split[0], k -> new TreeSet<>()).add(split[1]);
    }

    private void checkValidFilterEntry(String filter, String[] split) {
        String message = null;
        if (split.length != 2) {
            message = "Invalid read group filter: " + filter;
        } else if (split[0].length() != 2) {
            message = "Tag is not two characters: " + filter;
        }

        if (message != null) {
            message += ", format is <TAG>:<SUBSTRING>";
            throw new UserException(message);
        }
    }

    @Override
    public boolean test( final GATKRead read ) {
        final SAMReadGroupRecord readGroup = ReadUtils.getSAMReadGroupRecord(read, samHeader);
        if ( readGroup == null ) {
            return true;
        }

        for (final String attributeType : blacklistEntries.keySet()) {

            final String attribute;
            if (SAMReadGroupRecord.READ_GROUP_ID_TAG.equals(attributeType) || SAMTag.RG.name().equals(attributeType)) {
                attribute = readGroup.getId();
            } else {
                attribute = readGroup.getAttribute(attributeType);
            }
            if (attribute != null && blacklistEntries.get(attributeType).contains(attribute)) {
                return false;
            }
        }

        return true;
    }

}
