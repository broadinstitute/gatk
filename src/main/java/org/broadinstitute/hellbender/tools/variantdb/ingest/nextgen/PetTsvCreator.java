package org.broadinstitute.hellbender.tools.variantdb.ingest.nextgen;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.variantdb.ingest.IngestConstants;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

public final class PetTsvCreator {
    private static final Logger logger = LogManager.getLogger(PetTsvCreator.class);

    private SimpleXSVWriter petWriter = null;
    private final String sampleId;
    private SimpleInterval previousInterval;
    private final SAMSequenceDictionary seqDictionary;
    private final GQStateEnum gqStateToIgnore;
    private GenomeLocSortedSet coverageLocSortedSet;

    public PetTsvCreator(String sampleName, String sampleId, Path sampleDirectoryPath, SAMSequenceDictionary seqDictionary, GQStateEnum gqStateToIgnore) {
        this.sampleId = sampleId;
        this.seqDictionary = seqDictionary;
        this.gqStateToIgnore = gqStateToIgnore;
        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

        // If the pet directory inside it doesn't exist yet -- create it
        final String petDirectoryName = "pet";
        final Path petDirectoryPath = sampleDirectoryPath.resolve(petDirectoryName);
        final File petDirectory = new File(petDirectoryPath.toString());
        if (!petDirectory.exists()) {
            petDirectory.mkdir();
        }
       try {
            // Create a pet file to go into the pet dir for _this_ sample
            final String petOutputName = sampleName + petDirectoryName + IngestConstants.FILETYPE;
            final Path petOutputPath = petDirectoryPath.resolve(petOutputName);
            // write header to it
            List<String> petHeader = PetTsvCreator.getHeaders();
            petWriter = new SimpleXSVWriter(petOutputPath, IngestConstants.SEPARATOR);
            petWriter.setHeaderLine(petHeader);
        } catch (final IOException e) {
            throw new UserException("Could not create pet outputs", e);
        }

    }

    /**
     * Expected headers for the Position Table (PET)
     */
    public enum PetFieldEnum {
        position,
        sample,
        state,
    }

    public enum GQStateEnum {
        VARIANT("v"),
        STAR("*"),
        ZERO("0"),
        TEN("1"),
        TWENTY("2"),
        THIRTY("3"),
        FORTY("4"),
        FIFTY("5"),
        SIXTY("6"),
        MISSING("m"),
        UNKNOWN("u");

        String value;

        GQStateEnum(String value) {
            this.value = value;
        }

    }

    public void apply(VariantContext variant, List<GenomeLoc> intervalsToWrite) {
        boolean firstInterval = true;
        final String variantChr = variant.getContig();

        for (GenomeLoc genomeLoc : intervalsToWrite) {

            int start = Math.max(genomeLoc.getStart(), variant.getStart());
            int end = Math.min(genomeLoc.getEnd(), variant.getEnd());

            // TODO throw an error if start and end are the same?

            // for each of the reference blocks with the GQ to discard, keep track of the positions for the missing insertions
            if (PetTsvCreator.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(gqStateToIgnore)) {
                // add interval to "covered" intervals
                setCoveredInterval(variantChr, start, end);
            }

            // create PET output if the reference block's GQ is not the one to discard or its a variant
            if (!variant.isReferenceBlock() || !PetTsvCreator.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(gqStateToIgnore)) {

                // add interval to "covered" intervals
                setCoveredInterval(variantChr, start, end);

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
                    petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
                }
            }
            firstInterval = false;
        }

    }

    public void writeMissingIntervals(GenomeLocSortedSet intervalArgumentGenomeLocSortedSet) {
        GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(coverageLocSortedSet);
        logger.info("MISSING_GREP_HERE:" + uncoveredIntervals.coveredSize());
        logger.info("MISSING_PERCENTAGE_GREP_HERE:" + (1.0 * uncoveredIntervals.coveredSize()) / intervalArgumentGenomeLocSortedSet.coveredSize());
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            final String contig = genomeLoc.getContig();
            // write the position to the XSV
            for (List<String> TSVLineToCreatePet : PetTsvCreator.createMissingTSV(
                    SchemaUtils.encodeLocation(contig, genomeLoc.getStart()),
                    SchemaUtils.encodeLocation(contig, genomeLoc.getEnd()),
                    sampleId
            )) {
                petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
            }
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

        if (!variant.isReferenceBlock()) {
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(start));
            row.add(sampleId);
            row.add(GQStateEnum.VARIANT.value);
            rows.add(row);

            //if variant is variant and has additional positions--must be a deletion: add `*` state
            for (long i = start + 1 ; i <= end; i++){
                row = new ArrayList<>();
                row.add(String.valueOf(i));
                row.add(sampleId);
                row.add(GQStateEnum.STAR.value);
                rows.add(row);
            }
        } else {
            // TODO check in the tool to make sure it's only one sample
            GQStateEnum state = getGQStateEnum(variant.getGenotype(0).getGQ());

            for (long position = start; position <= end; position++){ // break up ref blocks
                List<String> row = new ArrayList<>();

                row.add(String.valueOf(position));
                row.add(sampleId);
                row.add(state.value);
                rows.add(row);
            }
        }

        return rows;
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
            row.add(GQStateEnum.STAR.value);
            rows.add(row);
        }

        return rows;
    }

    public static List<List<String>> createMissingTSV(long start, long end, String sampleName) {
        List<List<String>> rows = new ArrayList<>();

        for (long position = start; position <= end; position ++){
            List<String> row = new ArrayList<>();
            row.add(String.valueOf(position));
            row.add(sampleName);
            row.add(GQStateEnum.MISSING.value);
            rows.add(row);
        }

        return rows;
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

    public static List<String> getHeaders() {
        return Arrays.stream(PetFieldEnum.values()).map(String::valueOf).collect(Collectors.toList());
    }

    public void closeTool() {
        if (petWriter != null) {
            try {
                petWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }

    }
}
