package org.broadinstitute.hellbender.tools.variantdb.ingest;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ExampleProgramGroup;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.*;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.tsv.SimpleXSVWriter;

import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.nio.file.Path;
import java.io.File;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class IngestVariantWalker extends VariantWalker {
    static final Logger logger = LogManager.getLogger(IngestVariantWalker.class);

    public static final char SEPARATOR = '\t';
    public static final String FILETYPE = ".tsv";
    private SimpleXSVWriter vetWriter = null;
    private SimpleXSVWriter petWriter = null;
    private SimpleXSVWriter sampleMetadataWriter = null;
    private IngestTSVCreator tsvCreator;

    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;
    private GenomeLocSortedSet coverageLocSortedSet;
    private SimpleInterval previousInterval;
    private String sampleName;
    private String sampleId;
    private List<SimpleInterval> userIntervals;

    // To determine which directory (and ultimately table) the sample's data will go into
    // Since tables have a limited number of samples (default is 4k)
    public int getSampleDirectoryNumber(String sampleId, int sampleMod) { // this is based on sample id
        // sample ids 1-4000 will go in directory 001
        int sampleIdInt = Integer.valueOf(sampleId); // TODO--should sampleId just get refactored as a long?
        // subtract 1 from the sample id to make it 1-index (or do we want to 0-index?) and add 1 to the dir
        int directoryNumber = Math.floorDiv((sampleIdInt - 1), sampleMod) + 1; // TODO omg write some unit tests
        return directoryNumber;
    }

    // Inside the parent directory, a directory for each chromosome will be created, with a pet directory and vet directory in each one.
    // Each pet and vet directory will hold all of the pet and vet tsvs for each sample
    // A metadata directory will be created, with a metadata tsv for each sample

    @Argument(fullName = "vet-pet-out-path",
            shortName = "VPO",
            doc = "Path to the directory where the variants TSVs and positions expanded TSVs should be written")
    public GATKPathSpecifier parentOutputDirectory = null;
    public Path parentDirectory = null;

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public IngestPetCreation.GQStateEnum gqStateToIgnore = IngestPetCreation.GQStateEnum.SIXTY;

    @Argument(fullName = "sample-name-mapping",
            shortName = "SNM",
            doc = "Sample name to sample id mapping",
            optional = true)
    public File sampleMap;

    @Argument(fullName = "is-array",
            shortName = "IA",
            doc = "Flag if the input vcf is an array",
            optional = true)
    public Boolean isArray = false;

    @Argument(fullName = "mode",
            shortName = "MO",
            doc = "Type of sample. Default is EXOMES. Valid options are ARRAYS, EXOMES, GENOMES",
            optional = true)
    public CommonCode.ModeEnum mode = CommonCode.ModeEnum.EXOMES;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    private String refVersion = "37";

    @Argument(
            fullName = "NextId",
            shortName = "ID",
            doc = "Id to start numbering the samples with",
            optional = true)
    private Integer nextId;


//    @Override
//    public boolean requiresIntervals() {
//        return true; // TODO -- do I need to check the boolean flag on this?
//    }

    @Override
    public void onTraversalStart() {

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);

        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        sampleName = IngestUtils.getSampleName(inputVCFHeader);
        if (sampleMap == null) {
            sampleId = nextId.toString();
            nextId = nextId++;
        } else {
            sampleId = IngestUtils.getSampleId(sampleName, sampleMap);
        }

        // Mod the sample directories
        final int sampleMod = 4000; // TODO hardcoded for now--potentially an input param?
        int sampleDirectoryNumber = getSampleDirectoryNumber(sampleId, sampleMod);

        parentDirectory = parentOutputDirectory.toPath(); // TODO do we need this? More efficient way to do this?
        // If this sample set directory doesn't exist yet -- create it
        final String sampleDirectoryName = String.valueOf(sampleDirectoryNumber);
        final Path sampleDirectoryPath = parentDirectory.resolve(sampleDirectoryName);
        final File sampleDirectory = new File(sampleDirectoryPath.toString());
        if (! sampleDirectory.exists()){
            sampleDirectory.mkdir();
        }

        // If the metadata directory inside it doesn't exist yet, create it
        final String metadataDirectoryName = "metadata";
        final Path metadataDirectoryPath = sampleDirectoryPath.resolve(metadataDirectoryName);
        final File sampleMetadataOutputDirectory = new File(metadataDirectoryPath.toString());
        if (! sampleMetadataOutputDirectory.exists()){
            sampleMetadataOutputDirectory.mkdir();
        }


        if (isArray) {
            tsvCreator = new IngestRawArrayTSVCreator(sampleName, sampleId, sampleDirectoryPath);
        } else {
            final SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();
            userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

            // To set up the missing positions
            final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
            intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));
            coverageLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));


            // If the pet directory inside it doesn't exist yet -- create it
            final String petDirectoryName = "pet";
            final Path petDirectoryPath = sampleDirectoryPath.resolve(petDirectoryName);
            final File petDirectory = new File(petDirectoryPath.toString());
            if (!petDirectory.exists()) {
                petDirectory.mkdir();
            }
            // If the vet directory inside it doesn't exist yet -- create it
            final String vetDirectoryName = "vet";
            final Path vetDirectoryPath = sampleDirectoryPath.resolve(vetDirectoryName);
            final File vetDirectory = new File(vetDirectoryPath.toString());
            if (!vetDirectory.exists()) {
                vetDirectory.mkdir();
            }

            // if the pet & vet tsvs don't exist yet -- create them
            try {
                // Create a pet file to go into the pet dir for _this_ sample
                final String petOutputName = sampleName + petDirectoryName + FILETYPE;
                final Path petOutputPath = petDirectoryPath.resolve(petOutputName);
                // write header to it
                List<String> petHeader = IngestPetCreation.getHeaders();
                petWriter = new SimpleXSVWriter(petOutputPath, SEPARATOR);
                petWriter.setHeaderLine(petHeader);
            } catch (final IOException e) {
                throw new UserException("Could not create pet outputs", e);
            }

            try {
                // Create a vet file to go into the pet dir for _this_ sample
                final String vetOutputName = sampleName + vetDirectoryName + FILETYPE;
                final Path vetOutputPath = vetDirectoryPath.resolve(vetOutputName);
                // write header to it
                List<String> vetHeader = IngestVetExomeCreation.getHeaders();
                vetWriter = new SimpleXSVWriter(vetOutputPath, SEPARATOR);
                vetWriter.setHeaderLine(vetHeader);
            } catch (final IOException e) {
                throw new UserException("Could not create vet outputs", e);
            }


        }


        // if the metadata tsvs don't exist yet -- create them
        try {
            // Create a metadata file to go into the metadata dir for _this_ sample
            // TODO--this should just be one file per sample set?
            final String sampleMetadataName = sampleName + metadataDirectoryName + FILETYPE;
            final Path sampleMetadataOutputPath = metadataDirectoryPath.resolve(sampleMetadataName);
            // write header to it
            List<String> sampleListHeader = IngestSampleListCreation.getHeaders();
            sampleMetadataWriter = new SimpleXSVWriter(sampleMetadataOutputPath, SEPARATOR);
            sampleMetadataWriter.setHeaderLine(sampleListHeader);
            String intervalListMd5 = "NA";

            if (!isArray) {
                // write values
                List<String> intervalList = userIntervals.stream().map(interval -> interval.toString())
                        .collect(Collectors.toList());
                String intervalListBlob = StringUtils.join(intervalList, ", ");
                intervalListMd5 = Utils.calcMD5(intervalListBlob);
            }
            final List<String> TSVLineToCreateSampleMetadata = IngestSampleListCreation.createSampleListRow(
                    sampleName,
                    sampleId,
                    intervalListMd5,
                    gqStateToIgnore);
            sampleMetadataWriter.getNewLineBuilder().setRow(TSVLineToCreateSampleMetadata).write();

        } catch (final IOException e) {
            throw new UserException("Could not create sample metadata outputs", e);
        }
    }

    public static void setCoveredInterval(IngestVariantWalker ingestVariantWalker, String variantChr, int start, int end) {
        // add interval to "covered" intervals
        // GenomeLocSortedSet will automatically merge intervals that are overlapping when setting `mergeIfIntervalOverlaps`
        // to true.  In a GVCF most blocks are adjacent to each other so they wouldn't normally get merged.  We check
        // if the current record is adjacent to the previous record and "overlap" them if they are so our set is as
        // small as possible while still containing the same bases.
        final SimpleInterval variantInterval = new SimpleInterval(variantChr, start, end);
        final int intervalStart = (ingestVariantWalker.previousInterval != null && ingestVariantWalker.previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                ingestVariantWalker.previousInterval.getStart() : variantInterval.getStart();
        final int intervalEnd = (ingestVariantWalker.previousInterval != null && ingestVariantWalker.previousInterval.overlapsWithMargin(variantInterval, 1)) ?
                Math.max(ingestVariantWalker.previousInterval.getEnd(), variantInterval.getEnd()) : variantInterval.getEnd();

        final GenomeLoc possiblyMergedGenomeLoc = ingestVariantWalker.coverageLocSortedSet.getGenomeLocParser().createGenomeLoc(variantInterval.getContig(), intervalStart, intervalEnd);
        ingestVariantWalker.coverageLocSortedSet.add(possiblyMergedGenomeLoc, true);
        ingestVariantWalker.previousInterval = new SimpleInterval(possiblyMergedGenomeLoc);
    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

        if (isArray) {
            tsvCreator.apply(variant, readsContext, referenceContext, featureContext);
            return;
        }

        // get the intervals this variant covers
        final GenomeLoc variantGenomeLoc = intervalArgumentGenomeLocSortedSet.getGenomeLocParser().createGenomeLoc(variant.getContig(), variant.getStart(), variant.getEnd());
        final List<GenomeLoc> intervalsToWrite = intervalArgumentGenomeLocSortedSet.getOverlapping(variantGenomeLoc);

        if (intervalsToWrite.size() == 0){
            throw new IllegalStateException("There are no intervals being covered by this variant, something went wrong with interval parsing");
        }

        // take the first interval(assuming this is returned in order) and make sure if its a variant, that it starts at/after the interval start
        // we are going to ignore any deletions that start before an interval.
        if (!variant.isReferenceBlock() && intervalsToWrite.get(0).getStart() > variant.getStart()){
            return;
        }

        // if the only alt allele for a variant is `*`, we ignore it
        if (!variant.isReferenceBlock() &&  variant.getAlternateAlleles().size() == 2 && variant.hasAlternateAllele(Allele.SPAN_DEL)){
            return;
        }

        final String variantChr = variant.getContig();

        // create VET output
        if (!variant.isReferenceBlock()) {
            int start = variant.getStart();
            int end = variant.getEnd();
            // check to see if this is an array
            // else, it must be an exome or genome!
            final List<String> TSVLineToCreateVet = IngestVetExomeCreation.createVariantRow(
                    SchemaUtils.encodeLocation(variantChr, start),
                    variant,
                    sampleId
            );

            // write the variant to the XSV
            SimpleXSVWriter.LineBuilder vetLine = vetWriter.getNewLineBuilder();
            vetLine.setRow(TSVLineToCreateVet);
            vetLine.write();
        }
        boolean firstInterval = true;
        for (GenomeLoc genomeLoc : intervalsToWrite) {

            int start = Math.max(genomeLoc.getStart(), variant.getStart());
            int end = Math.min(genomeLoc.getEnd(), variant.getEnd());

            // TODO throw an error if start and end are the same?

            // for each of the reference blocks with the GQ to discard, keep track of the positions for the missing insertions
            if (IngestPetCreation.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(gqStateToIgnore)) {
                // add interval to "covered" intervals
                setCoveredInterval(this, variantChr, start, end);
            }

            // create PET output if the reference block's GQ is not the one to discard or its a variant
            if (!variant.isReferenceBlock() || !IngestPetCreation.getGQStateEnum(variant.getGenotype(0).getGQ()).equals(gqStateToIgnore)) {

                // add interval to "covered" intervals
                setCoveredInterval(this, variantChr, start, end);

                List<List<String>> TSVLinesToCreatePet;
                // handle deletions that span across multiple intervals
                if (!firstInterval && !variant.isReferenceBlock()) {
                    TSVLinesToCreatePet = IngestPetCreation.createSpanDelRows(
                            SchemaUtils.encodeLocation(variantChr, start),
                            SchemaUtils.encodeLocation(variantChr, end),
                            variant,
                            sampleId
                    );
                } else {
                    TSVLinesToCreatePet = IngestPetCreation.createPositionRows(
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

    @Override
    public Object onTraversalSuccess() {
        if (isArray) { return 0; }
        final GenomeLocSortedSet uncoveredIntervals = intervalArgumentGenomeLocSortedSet.subtractRegions(coverageLocSortedSet);
        logger.info("MISSING_GREP_HERE:" + uncoveredIntervals.coveredSize());
        logger.info("MISSING_PERCENTAGE_GREP_HERE:" + (1.0 * uncoveredIntervals.coveredSize()) / intervalArgumentGenomeLocSortedSet.coveredSize());
        for (GenomeLoc genomeLoc : uncoveredIntervals) {
            final String contig = genomeLoc.getContig();
            // write the position to the XSV
            for (List<String> TSVLineToCreatePet : IngestPetCreation.createMissingTSV(
                    SchemaUtils.encodeLocation(contig, genomeLoc.getStart()),
                    SchemaUtils.encodeLocation(contig, genomeLoc.getEnd()),
                    sampleId
            )) {
                petWriter.getNewLineBuilder().setRow(TSVLineToCreatePet).write();
            }
        }
        return 0;
    }

    @Override
    public void closeTool() {
        if (isArray) {
            tsvCreator.closeTool();
        }
        if (vetWriter != null) {
            try {
                vetWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close VET writer", e);
            }
        }
        if (petWriter != null) {
            try {
                petWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close PET writer", e);
            }
        }
        if (sampleMetadataWriter != null) {
            try {
                sampleMetadataWriter.close();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close Sample Metadata writer", e);
            }
        }
    }
}
