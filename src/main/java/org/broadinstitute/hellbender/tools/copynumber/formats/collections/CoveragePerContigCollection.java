package org.broadinstitute.hellbender.tools.copynumber.formats.collections;

import org.apache.commons.collections4.ListUtils;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.DetermineGermlineContigPloidy;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.Metadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.records.CoveragePerContig;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents a sequence dictionary and total coverage over each contig in an ordered set associated with a cohort of named samples.
 * Should only be used to write temporary files in {@link DetermineGermlineContigPloidy}.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CoveragePerContigCollection extends AbstractRecordCollection<LocatableMetadata, CoveragePerContig> {
    private static final String SAMPLE_NAME_TABLE_COLUMN = "SAMPLE_NAME";

    public CoveragePerContigCollection(final LocatableMetadata metadata,
                                       final List<CoveragePerContig> coveragePerContigs,
                                       final List<String> contigs) {
        super(
                metadata,
                coveragePerContigs,
                new TableColumnCollection(ListUtils.union(Collections.singletonList(SAMPLE_NAME_TABLE_COLUMN), contigs)),
                dataLine -> new CoveragePerContig(
                        dataLine.get(SAMPLE_NAME_TABLE_COLUMN),
                        contigs.stream().collect(Collectors.toMap(
                                Function.identity(),
                                dataLine::getInt,
                                (u, v) -> {
                                    throw new GATKException.ShouldNeverReachHereException("Cannot have duplicate contigs.");
                                },   //contigs should already be distinct
                                LinkedHashMap::new))),
                (coveragePerContig, dataLine) -> {
                    dataLine.append(coveragePerContig.getSampleName());
                    contigs.stream().map(coveragePerContig::getCoverage).forEach(dataLine::append);
                });
    }

    @Override
    Metadata.Type getMetadataType() {
        return Metadata.Type.LOCATABLE;
    }
}
