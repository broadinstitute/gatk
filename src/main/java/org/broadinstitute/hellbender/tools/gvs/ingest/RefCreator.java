package org.broadinstitute.hellbender.tools.gvs.ingest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.gvs.common.CommonCode;
import org.broadinstitute.hellbender.tools.gvs.common.GQStateEnum;
import org.broadinstitute.hellbender.tools.gvs.common.IngestConstants;
import org.broadinstitute.hellbender.tools.gvs.common.SchemaUtils;
import org.broadinstitute.hellbender.utils.GenomeLoc;
import org.broadinstitute.hellbender.utils.GenomeLocParser;
import org.broadinstitute.hellbender.utils.GenomeLocSortedSet;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.BigQueryUtils;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public final class RefCreator {
    private static final Logger logger = LogManager.getLogger(RefCreator.class);

    private final CommonCode.OutputType outputType;

    private RefRangesWriter refRangesWriter = null;

    private final boolean writeReferenceRanges;
    private final String sampleId;
    private SimpleInterval previousInterval;
    private final Set<GQStateEnum> gqStatesToIgnore = new HashSet<>();
    private final GenomeLocSortedSet coverageLocSortedSet;
    private static final String PREFIX_SEPARATOR = "_";
    private final static String REF_RANGES_FILETYPE_PREFIX = "ref_ranges_";

    public static boolean doRowsExistFor(CommonCode.OutputType outputType, String projectId, String datasetName, String tableNumber, String sampleId) {
        if (outputType != CommonCode.OutputType.BQ) return false;
        return BigQueryUtils.doRowsExistFor(projectId, datasetName, REF_RANGES_FILETYPE_PREFIX + tableNumber, sampleId);
    }

    public RefCreator(String sampleIdentifierForOutputFileName, String sampleId, String tableNumber, SAMSequenceDictionary seqDictionary, GQStateEnum gqStateToIgnore, final boolean dropAboveGqThreshold, final File outputDirectory, final CommonCode.OutputType outputType, final boolean writeReferenceRanges, final String projectId, final String datasetName) {
        this.sampleId = sampleId;
        this.outputType = outputType;
        this.writeReferenceRanges = writeReferenceRanges;

        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

        try {
            if (writeReferenceRanges) {
                final File refOutputFile = new File(outputDirectory, REF_RANGES_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + "." + outputType.toString().toLowerCase());
                switch (outputType) {
                    case BQ:
                        if (projectId == null || datasetName == null) {
                            throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                        }
                        refRangesWriter = new RefRangesBQWriter(projectId, datasetName,REF_RANGES_FILETYPE_PREFIX + tableNumber);
                        break;
                    case TSV:
                        refRangesWriter = new RefRangesTsvWriter(refOutputFile.getCanonicalPath());
                        break;
                    case AVRO:
                        refRangesWriter = new RefRangesAvroWriter(refOutputFile.getCanonicalPath());
                        break;
                }
            }
        } catch (final IOException ioex) {
            throw new UserException("Could not create reference range outputs", ioex);
        }

        this.gqStatesToIgnore.add(gqStateToIgnore);
        if (dropAboveGqThreshold) {
            this.gqStatesToIgnore.addAll(getGQStateEnumGreaterThan(gqStateToIgnore));
        }
    }

    public void apply(VariantContext variant, List<GenomeLoc> intervalsToWrite) throws IOException {
        final String variantChr = variant.getContig();

        for (GenomeLoc genomeLoc : intervalsToWrite) {

            int start = Math.max(genomeLoc.getStart(), variant.getStart());
            int end = Math.min(genomeLoc.getEnd(), variant.getEnd());

            // TODO throw an error if start and end are the same?

            // for each of the reference blocks with the GQ to discard, keep track of the positions for the missing insertions
            if (this.gqStatesToIgnore.contains(getGQStateEnum(variant.getGenotype(0).getGQ()))) {
                // add interval to "covered" intervals
                setCoveredInterval(variantChr, start, end);
            }

            // Create output if the reference block's GQ is not the one to discard, or it's a variant.
            if (!variant.isReferenceBlock() || !this.gqStatesToIgnore.contains(RefCreator.getGQStateEnum(variant.getGenotype(0).getGQ()))) {

                // add interval to "covered" intervals
                setCoveredInterval(variantChr, start, end);

                // if we are writing ref ranges, and this is a reference block, write it!
                if (writeReferenceRanges) {
                    if (variant.isReferenceBlock()) {
                        // break up reference blocks to be no longer than MAX_REFERENCE_BLOCK_SIZE
                        int localStart = start;
                        while ( localStart <= end ) {
                            int length = Math.min(end - localStart + 1, IngestConstants.MAX_REFERENCE_BLOCK_BASES);
                            refRangesWriter.write(SchemaUtils.encodeLocation(variantChr, localStart),
                                    Long.parseLong(sampleId),
                                    length,
                                    getGQStateEnum(variant.getGenotype(0).getGQ()).getValue()
                            );
                            localStart = localStart + length ;
                        }

                    // Write out no-calls as a single-base GQ0 reference.
                    } else if (CreateVariantIngestFiles.isNoCall(variant)) {
                        refRangesWriter.write(SchemaUtils.encodeLocation(variantChr, start),
                                Long.parseLong(sampleId),
                                1,
                                GQStateEnum.ZERO.getValue()
                        );
                    }
                }
            }
        }
    }

    public void writeMissingIntervals(GenomeLocSortedSet intervalArgumentGenomeLocSortedSet) throws IOException {
        GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(coverageLocSortedSet);
        logger.info("MISSING_GREP_HERE:" + uncoveredIntervals.coveredSize());
        logger.info("MISSING_PERCENTAGE_GREP_HERE:" + (1.0 * uncoveredIntervals.coveredSize()) / intervalArgumentGenomeLocSortedSet.coveredSize());
        // for each block of uncovered locations
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            final String contig = genomeLoc.getContig();
            // write all positions in this block
            writeMissingPositions(
                    SchemaUtils.encodeLocation(contig, genomeLoc.getStart()),
                    SchemaUtils.encodeLocation(contig, genomeLoc.getEnd()));
        }
    }

    private void setCoveredInterval(String variantChr, int start, int end) {
        // add interval to "covered" intervals
        // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
        // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
        // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
        // small as possible while still containing the same bases.
        final SimpleInterval variantInterval = new SimpleInterval(variantChr, start, end);

        boolean overlapping = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1));
        final int intervalStart = overlapping ? previousInterval.getStart() : variantInterval.getStart();
        final int intervalEnd = overlapping ? Math.max(previousInterval.getEnd(), variantInterval.getEnd()) : variantInterval.getEnd();

        final GenomeLoc possiblyMergedGenomeLoc = coverageLocSortedSet.getGenomeLocParser().createGenomeLoc(variantInterval.getContig(), intervalStart, intervalEnd);
        coverageLocSortedSet.add(possiblyMergedGenomeLoc, true);
        previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
    }

    public void writeMissingPositions(long start, long end) throws IOException {
        if (writeReferenceRanges) {
            // break up missing blocks to be no longer than MAX_REFERENCE_BLOCK_SIZE
            long localStart = start;
            while ( localStart <= end ) {
                int length = (int) Math.min(end - localStart + 1, IngestConstants.MAX_REFERENCE_BLOCK_BASES);
                refRangesWriter.write(localStart,
                        Long.parseLong(sampleId),
                        length,
                        GQStateEnum.MISSING.getValue()
                );
                localStart = localStart + length ;
            }
        }
    }

    public static GQStateEnum getGQStateEnum(int GQ){
        if (GQ < 10) {
            return GQStateEnum.ZERO;
        } else if (GQ < 20) {
            return GQStateEnum.TEN;
        } else if (GQ < 30) {
            return GQStateEnum.TWENTY;
        } else if (GQ < 40) {
            return GQStateEnum.THIRTY;
        } else if (GQ < 50) {
            return GQStateEnum.FORTY;
        } else if (GQ < 60) {
            return GQStateEnum.FIFTY;
        } else {
            return GQStateEnum.SIXTY;
        }
    }

    // this is ugly.... I think we need to rework the enum to better handle the new use cases
    // but just getting this going.
    public static Set<GQStateEnum> getGQStateEnumGreaterThan(GQStateEnum s){
        Set<GQStateEnum> ret = new HashSet<>();

        switch (s) {
            case ZERO:
                ret.add(GQStateEnum.TEN);
                ret.add(GQStateEnum.TWENTY);
                ret.add(GQStateEnum.THIRTY);
                ret.add(GQStateEnum.FORTY);
                ret.add(GQStateEnum.FIFTY);
                ret.add(GQStateEnum.SIXTY);
                break;
            case TEN:
                ret.add(GQStateEnum.TWENTY);
                ret.add(GQStateEnum.THIRTY);
                ret.add(GQStateEnum.FORTY);
                ret.add(GQStateEnum.FIFTY);
                ret.add(GQStateEnum.SIXTY);
                break;
            case TWENTY:
                ret.add(GQStateEnum.THIRTY);
                ret.add(GQStateEnum.FORTY);
                ret.add(GQStateEnum.FIFTY);
                ret.add(GQStateEnum.SIXTY);
                break;
            case THIRTY:
                ret.add(GQStateEnum.FORTY);
                ret.add(GQStateEnum.FIFTY);
                ret.add(GQStateEnum.SIXTY);
                break;
            case FORTY:
                ret.add(GQStateEnum.FIFTY);
                ret.add(GQStateEnum.SIXTY);
                break;
            case FIFTY:
                ret.add(GQStateEnum.SIXTY);
                break;
        }

        return ret;
    }

    public void commitData() {
        if (outputType == CommonCode.OutputType.BQ) {
            if (writeReferenceRanges && refRangesWriter != null) {
                refRangesWriter.commitData();
            }
        }
    }

    public void closeTool() {
        try {
            if (refRangesWriter != null) refRangesWriter.close();
        } catch (final Exception e) {
            throw new IllegalArgumentException("Couldn't close reference ranges writer", e);
        }
    }
}
