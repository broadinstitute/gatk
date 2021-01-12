package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.*;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.IngestUtils;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.io.File;
import java.io.IOException;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Exome and Genome Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateVariantIngestFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateVariantIngestFiles.class);

    private PetTsvCreator petTsvCreator;
    private VetTsvCreator vetTsvCreator;
    private MetadataTsvCreator metadataTsvCreator;
    private GenomeLocSortedSet intervalArgumentGenomeLocSortedSet;

    private String sampleName;
    private String sampleId;
    private List<SimpleInterval> userIntervals;

    // Inside the parent directory, a directory for each chromosome will be created, with a pet directory and vet directory in each one.
    // Each pet and vet directory will hold all of the pet and vet tsvs for each sample
    // A metadata directory will be created, with a metadata tsv for each sample

//    @Argument(fullName = "output-path",
//            shortName = "VPO",
//            doc = "Path to the directory where the variants TSVs and positions expanded TSVs should be written")
//    public GATKPathSpecifier parentOutputDirectory = null;
//    public Path parentDirectory = null;

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public PetTsvCreator.GQStateEnum gqStateToIgnore = null;

    @Argument(fullName = "ignore-above-gq-threshold",
    shortName = "GTIG",
    doc = "in addition to dropping the gq block specified by ref-block-gq-to-ignore, also drop higher gq blocks",
    optional = true)
    public boolean dropAboveGqThreshold = false;

    @Argument(fullName = "emit-spanning-deletions-in-vet",
    shortName = "esdav",
    doc = "should spanning deletions be emitted as variants?",
    optional = true)
    public boolean emitSpanningDeletionsInVet = true;

    @Argument(fullName = "sample-name-mapping",
            shortName = "SNM",
            doc = "Sample name to sample id mapping",
            optional = true)
    public File sampleMap;

    @Argument(fullName = "sample-id",
            shortName = "SID",
            doc = "Sample id",
            optional = true)
    public Long sampleIdParam;

    @Argument(fullName = "mode",
            shortName = "MO",
            doc = "Type of sample. Default is EXOMES. Valid options are EXOMES, GENOMES",
            optional = true)
    public CommonCode.ModeEnum mode = CommonCode.ModeEnum.EXOMES;

    @Argument(fullName = "output-type", 
            shortName = "ot", 
            doc = "[Experimental] Output file format: TSV, ORC or PARQUET [default=TSV].", 
            optional = true)
    public CommonCode.OutputType outputType = CommonCode.OutputType.TSV;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    private String refVersion = "37";

    @Argument(
            fullName = "output-directory",
            doc = "directory for output tsv files",
            optional = true)
    private File outputDir = new File(".");


    @Override
    public boolean requiresIntervals() {
        return true; // TODO -- do I need to check the boolean flag on this?
    }

    @Override
    public void onTraversalStart() {
        //set up output directory
        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);

        // TODO should we reuse the SampleList class or move these methods there?
        // TODO if you change here, also change in CreateArrayIngestFiles
        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        sampleName = IngestUtils.getSampleName(inputVCFHeader);
        if (sampleIdParam == null && sampleMap == null) {
            throw new IllegalArgumentException("One of sample-id or sample-name-mapping must be specified");
        }
        if (sampleIdParam != null) {
            sampleId = String.valueOf(sampleIdParam);
        } else {
            sampleId = IngestUtils.getSampleId(sampleName, sampleMap);
        }


        // Mod the sample directories
        int sampleTableNumber = IngestUtils.getTableNumber(sampleId, IngestConstants.partitionPerTable);
        String tableNumberPrefix = String.format("%03d_", sampleTableNumber);

//        parentDirectory = parentOutputDirectory.toPath(); // TODO do we need this? More efficient way to do this?
//        final Path sampleDirectoryPath = IngestUtils.createSampleDirectory(parentDirectory, sampleDirectoryNumber);
        metadataTsvCreator = new MetadataTsvCreator(sampleName, sampleId, tableNumberPrefix, outputDir);
        metadataTsvCreator.createRow(sampleName, sampleId, userIntervals, gqStateToIgnore);

        // To set up the missing positions
        SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();
        userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));

        petTsvCreator = new PetTsvCreator(sampleName, sampleId, tableNumberPrefix, seqDictionary, gqStateToIgnore, dropAboveGqThreshold, outputDir, outputType);
        switch (mode) {
            case EXOMES:
            case GENOMES:
                vetTsvCreator = new VetTsvCreator(sampleName, sampleId, tableNumberPrefix, outputDir);
                break;
            case ARRAYS:
                throw new UserException.BadInput("To ingest Array data, use CreateArrayIngestFiles tool.");
        }


    }

    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {

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
        boolean isSpanDelVariant = !variant.isReferenceBlock() &&  variant.getAlternateAlleles().size() == 2 && variant.hasAlternateAllele(Allele.SPAN_DEL);

        // create VET output
        if (!variant.isReferenceBlock()) {
            // write VET entries for spanning deletion entires only if emitSpanningDeletionsInVet is enabled
            if (!isSpanDelVariant || emitSpanningDeletionsInVet) {
                vetTsvCreator.apply(variant, readsContext, referenceContext, featureContext);
            }
        }

        try {
            petTsvCreator.apply(variant, intervalsToWrite, isSpanDelVariant);
        } catch (IOException ioe) {
            throw new GATKException("Error writing PET", ioe);
        }

    }

    @Override
    public Object onTraversalSuccess() {
        try {
            petTsvCreator.writeMissingIntervals(intervalArgumentGenomeLocSortedSet);
        } catch (IOException ioe) {
            throw new GATKException("Error writing missing intervals", ioe);
        }
        return 0;
    }

    @Override
    public void closeTool() {
        if (petTsvCreator != null) {
            petTsvCreator.closeTool();
        }
        if (vetTsvCreator != null) {
            vetTsvCreator.closeTool();;
        }
        if (metadataTsvCreator != null) {
            try {
                metadataTsvCreator.closeTool();
            } catch (final Exception e) {
                throw new IllegalArgumentException("Couldn't close Sample Metadata writer", e);
            }
        }
    }

}
