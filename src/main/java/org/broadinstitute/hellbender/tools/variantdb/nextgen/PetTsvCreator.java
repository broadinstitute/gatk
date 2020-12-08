package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public final class PetTsvCreator {
    private static final Logger logger = LogManager.getLogger(PetTsvCreator.class);

    private CommonCode.OutputType outputType;

    private SimpleXSVWriter petTsvWriter = null;
    private PetTsvWriter petTsv2Writer = null;
    private PetOrcWriter petOrcWriter = null;
    private PetAvroWriter petAvroWriter = null;
    private PetParquetWriter petParquetWriter = null;

    private final String sampleId;
    private SimpleInterval previousInterval;
    private final SAMSequenceDictionary seqDictionary;
    private final Set<GQStateEnum> gqStatesToIgnore = new HashSet<GQStateEnum>();
    private GenomeLocSortedSet coverageLocSortedSet;
    private static String PET_FILETYPE_PREFIX = "pet_";

    public PetTsvCreator(String sampleName, String sampleId, String tableNumberPrefix, SAMSequenceDictionary seqDictionary, GQStateEnum gqStateToIgnore, final boolean dropAboveGqThreshold, final File outputDirectory, final CommonCode.OutputType outputType) {
        this.sampleId = sampleId;
        this.seqDictionary = seqDictionary;
        this.outputType = outputType;

        coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));

       try {
            final File petOutputFile = new File(outputDirectory, PET_FILETYPE_PREFIX + tableNumberPrefix + sampleName  + "." + outputType.toString().toLowerCase());
            switch (outputType) {
                case TSV:                
                    List<String> petHeader = PetTsvCreator.getHeaders();
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
        } catch (final IOException e) {
            throw new UserException("Could not create pet outputs", e);
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

        public String getValue() {
            return value;
        }

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
            if (!variant.isReferenceBlock() || !this.gqStatesToIgnore.contains(PetTsvCreator.getGQStateEnum(variant.getGenotype(0).getGQ()))) {

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
                    long location = Long.parseLong(TSVLineToCreatePet.get(0));
                    long sampleId = Long.parseLong(TSVLineToCreatePet.get(1));
                    String state = TSVLineToCreatePet.get(2);

                    switch (outputType) {
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
            firstInterval = false;
        }

    }

    public void writeMissingIntervals(GenomeLocSortedSet intervalArgumentGenomeLocSortedSet) throws IOException {
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
                long location = Long.parseLong(TSVLineToCreatePet.get(0));
                long sampleId = Long.parseLong(TSVLineToCreatePet.get(1));
                String state = TSVLineToCreatePet.get(2);

                switch (outputType) {
                    case TSV:                
                        petTsvWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
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

    public void closeTool() {
        try {
            switch (outputType) {
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
        } catch (final Exception e) {
            throw new IllegalArgumentException("Couldn't close PET writer", e);
        }
    }
}
