package org.broadinstitute.hellbender.tools.gvs.ingest;

import com.google.protobuf.Descriptors;
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
import org.broadinstitute.hellbender.utils.bigquery.PendingBQWriter;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;
import org.json.JSONObject;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;
import java.util.stream.Collectors;


public final class RefCreator {
    private static final Logger logger = LogManager.getLogger(RefCreator.class);

    private CommonCode.OutputType outputType;

    private SimpleXSVWriter petTsvWriter = null;
    private PetTsvWriter petTsv2Writer = null;
    private PetOrcWriter petOrcWriter = null;
    private PetAvroWriter petAvroWriter = null;
    private PetParquetWriter petParquetWriter = null;
    private PendingBQWriter petBQJsonWriter = null;

    private RefRangesWriter refRangesWriter = null;

    private boolean writePetData;
    private boolean writeReferenceRanges;
    private final String sampleId;
    private SimpleInterval previousInterval;
    private final SAMSequenceDictionary seqDictionary;
    private final Set<GQStateEnum> gqStatesToIgnore = new HashSet<GQStateEnum>();
    private GenomeLocSortedSet coverageLocSortedSet;
    private static final String PET_FILETYPE_PREFIX = "pet_";
    private static final String PREFIX_SEPARATOR = "_";
    private final static String REF_RANGES_FILETYPE_PREFIX = "ref_ranges_";




    public RefCreator(String sampleIdentifierForOutputFileName, String sampleId, String tableNumber, SAMSequenceDictionary seqDictionary, GQStateEnum gqStateToIgnore, final boolean dropAboveGqThreshold, final File outputDirectory, final CommonCode.OutputType outputType, final boolean writePetData, final boolean writeReferenceRanges, final String projectId, final String datasetName) {
        this.sampleId = sampleId;
        this.seqDictionary = seqDictionary;
        this.outputType = outputType;
        this.writePetData = writePetData;
        this.writeReferenceRanges = writeReferenceRanges;

        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

        try {
            final File petOutputFile = new File(outputDirectory, PET_FILETYPE_PREFIX + tableNumber + PREFIX_SEPARATOR + sampleIdentifierForOutputFileName + "." + outputType.toString().toLowerCase());
            switch (outputType) {
                case BQ:
                    if (projectId == null || datasetName == null) {
                        throw new UserException("Must specify project-id and dataset-name when using BQ output mode.");
                    }
                    petBQJsonWriter = new PendingBQWriter(projectId, datasetName,PET_FILETYPE_PREFIX + tableNumber);
                    break;
                case TSV:
                    List<String> petHeader = RefCreator.getHeaders();
                    petTsvWriter = new SimpleXSVWriter(petOutputFile.toPath(), IngestConstants.SEPARATOR);
                    petTsvWriter.setHeaderLine(petHeader);
                    break;
                case TSV2:
                    petTsv2Writer = new PetTsvWriter(petOutputFile.getCanonicalPath());
                    break;
                case ORC:
                    petOrcWriter = new PetOrcWriter(petOutputFile.getCanonicalPath());
                    break;
                case AVRO:
                    petAvroWriter = new PetAvroWriter(petOutputFile.getCanonicalPath());
                    break;
                case PARQUET:
                    petParquetWriter = new PetParquetWriter(petOutputFile.getCanonicalPath());
            }

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
            throw new UserException("Could not create pet outputs", ioex);
        }

        this.gqStatesToIgnore.add(gqStateToIgnore);
        if (dropAboveGqThreshold) {
            this.gqStatesToIgnore.addAll(getGQStateEnumGreaterThan(gqStateToIgnore));
        }
    }

    /**
     * Expected headers for the Position Table (PET)
     */
    public enum PetFieldEnum {
        location,
        sample,
        state,
    }

    public void apply(VariantContext variant, List<GenomeLoc> intervalsToWrite) throws IOException {
        boolean firstInterval = true;
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

            // create PET output if the reference block's GQ is not the one to discard or its a variant
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

                    // write out no-calls as a single-base GQ0 reference (to match PET behavior)
                    } else if (CreateVariantIngestFiles.isNoCall(variant)) {
                        refRangesWriter.write(SchemaUtils.encodeLocation(variantChr, start),
                                Long.parseLong(sampleId),
                                1,
                                GQStateEnum.ZERO.getValue()
                        );
                    }
                }

                if (writePetData) {
                    List<List<String>> TSVLinesToCreatePet;
                    // handle deletions that span across multiple intervals
                    if (!firstInterval && !variant.isReferenceBlock()) {
                        TSVLinesToCreatePet = createSpanDelRows(
                                SchemaUtils.encodeLocation(variantChr, start),
                                SchemaUtils.encodeLocation(variantChr, end),
                                variant,
                                sampleId
                        );
                    } else {
                        TSVLinesToCreatePet = createRows(
                                SchemaUtils.encodeLocation(variantChr, start),
                                SchemaUtils.encodeLocation(variantChr, end),
                                variant,
                                sampleId
                        );

                    }


                    // write the position to the XSV
                    for (List<String> TSVLineToCreatePet : TSVLinesToCreatePet) {
                        long location = Long.parseLong(TSVLineToCreatePet.get(0));
                        long sampleId = Long.parseLong(TSVLineToCreatePet.get(1));
                        String state = TSVLineToCreatePet.get(2);

                        switch (outputType) {
                            case BQ:
                                try {
                                    petBQJsonWriter.addJsonRow(createJsonRow(location, sampleId, state));
                                } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                                    throw new IOException("BQ exception", ex);
                                }
                                break;
                            case TSV:
                                petTsvWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                                break;
                            case TSV2:
                                petTsv2Writer.addRow(location, sampleId, state);
                                break;
                            case ORC:
                                petOrcWriter.addRow(location, sampleId, state);
                                break;
                            case AVRO:
                                petAvroWriter.addRow(location, sampleId, state);
                                break;
                            case PARQUET:
                                petParquetWriter.addRow(location, sampleId, state);
                                break;
                        }
                    }
                }
            }
            firstInterval = false;
        }

    }

    public void writeMissingIntervals(GenomeLocSortedSet intervalArgumentGenomeLocSortedSet) throws IOException {
        GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(coverageLocSortedSet);
        logger.info("MISSING_GREP_HERE:" + uncoveredIntervals.coveredSize());
        logger.info("MISSING_PERCENTAGE_GREP_HERE:" + (1.0 * uncoveredIntervals.coveredSize()) / intervalArgumentGenomeLocSortedSet.coveredSize());
        // for each block of uncovered locations
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            final String contig = genomeLoc.getContig();
            // write all positions in this block to the pet output
            writeMissingPositions(
                    SchemaUtils.encodeLocation(contig, genomeLoc.getStart()),
                    SchemaUtils.encodeLocation(contig, genomeLoc.getEnd()),
                    sampleId
            );
        }
    }

    private void setCoveredInterval(String variantChr, int start, int end) {
        // add interval to "covered" intervals
        // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
        // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
        // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
        // small as possible while still containing the same bases.
        final SimpleInterval variantInterval = new SimpleInterval(variantChr, start, end);
        final int intervalStart = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                previousInterval.getStart() : variantInterval.getStart();
        final int intervalEnd = (previousInterval != null && previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                Math.max(previousInterval.getEnd(), variantInterval.getEnd()) : variantInterval.getEnd();

        final GenomeLoc possiblyMergedGenomeLoc = coverageLocSortedSet.getGenomeLocParser().createGenomeLoc(variantInterval.getContig(), intervalStart, intervalEnd);
        coverageLocSortedSet.add(possiblyMergedGenomeLoc, true);
        previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
    }

    public List<List<String>> createRows(final long start, final long end, final VariantContext variant, final String sampleId) {

        List<List<String>> rows = new ArrayList<>();

        // if the variant is no call, set the PET "state" to GQ ZERO
        if (CreateVariantIngestFiles.isNoCall(variant)) {
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(start));
            row.add(sampleId);
            row.add(GQStateEnum.ZERO.getValue());
            rows.add(row);
        } else if (!variant.isReferenceBlock()) {
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(start));
            row.add(sampleId);
            row.add(GQStateEnum.VARIANT.getValue());
            rows.add(row);

            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (long i = start + 1 ; i <= end; i++){
                row = new ArrayList<>();
                row.add(String.valueOf(i));
                row.add(sampleId);
                row.add(GQStateEnum.STAR.getValue());
                rows.add(row);
            }
        } else {
            // TODO check in the tool to make sure it's only one sample
            GQStateEnum state = getGQStateEnum(variant.getGenotype(0).getGQ());

            for (long position = start; position <= end; position++){ // break up ref blocks
                List<String> row = new ArrayList<>();

                row.add(String.valueOf(position));
                row.add(sampleId);
                row.add(state.getValue());
                rows.add(row);
            }
        }

        return rows;
    }

    private JSONObject createJsonRow(long location, long sampleId, String state) {
        JSONObject record = new JSONObject();
        record.put(SchemaUtils.LOCATION_FIELD_NAME, location);
        record.put(SchemaUtils.SAMPLE_ID_FIELD_NAME, sampleId);
        record.put(SchemaUtils.STATE_FIELD_NAME, state);
        return record;
    }


    public List<List<String>> createSpanDelRows(final long start, final long end, final VariantContext variant, final String sampleName) {
        if (variant.isReferenceBlock()){
            throw new IllegalStateException("Cannot create span deletion rows for a reference block");
        }

        List<List<String>> rows = new ArrayList<>();

        for (long position = start; position <= end; position++){ // break up ref blocks
            List<String> row = new ArrayList<>();

            row.add(String.valueOf(position));
            row.add(sampleName);
            row.add(GQStateEnum.STAR.getValue());
            rows.add(row);
        }

        return rows;
    }

    public void writeMissingPositions(long start, long end, String sampleName) throws IOException {
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

        if (writePetData) {
            writeMissingPetPositions(start, end, sampleName);
        }
    }

     public void writeMissingPetPositions(long start, long end, String sampleName) throws IOException {
        for (long position = start; position <= end; position++){
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(position));
            row.add(sampleName);
            row.add(GQStateEnum.MISSING.getValue());

            // TODO refactor - this only needs to be done for non-TSV outputTypes
            long location = Long.parseLong(row.get(0));
            long sampleId = Long.parseLong(row.get(1));
            String state = row.get(2);

            switch (outputType) {
                case TSV:
                    petTsvWriter.getNewLineBuilder().setRow(row).write();
                    break;
                case BQ:
                    try {
                        petBQJsonWriter.addJsonRow(createJsonRow(location, sampleId, state));
                    } catch (Descriptors.DescriptorValidationException | ExecutionException | InterruptedException ex) {
                        throw new IOException("BQ exception", ex);
                    }
                    break;
                case ORC:
                    petOrcWriter.addRow(location, sampleId, state);
                    break;
                case AVRO:
                    petAvroWriter.addRow(location, sampleId, state);
                    break;
                case PARQUET:
                    petParquetWriter.addRow(location, sampleId, state);
                    break;
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
        Set<GQStateEnum> ret = new HashSet<GQStateEnum>();

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

    public static List<String> getHeaders() {
        return Arrays.stream(PetFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void commitData() {
        if (outputType == CommonCode.OutputType.BQ) {
            if (writePetData && petBQJsonWriter != null) {
                petBQJsonWriter.flushBuffer();
                petBQJsonWriter.commitWriteStreams();
            }
            if (writeReferenceRanges && refRangesWriter != null) {
                refRangesWriter.commitData();
            }
        }
    }

    public void closeTool() {
        try {
            switch (outputType) {
                case BQ:
                    if (petBQJsonWriter != null) petBQJsonWriter.close();
                    break;
                case TSV:
                    if (petTsvWriter != null) petTsvWriter.close();
                    break;
                case ORC:
                    if (petOrcWriter != null) petOrcWriter.close();
                    break;
                case AVRO:
                    if (petAvroWriter != null) petAvroWriter.close();
                    break;
                case PARQUET:
                    if (petParquetWriter != null) petParquetWriter.close();
                    break;
            }

            if (refRangesWriter != null) refRangesWriter.close();

        } catch (final Exception e) {
            throw new IllegalArgumentException("Couldn't close PET writer", e);
        }
    }
}
