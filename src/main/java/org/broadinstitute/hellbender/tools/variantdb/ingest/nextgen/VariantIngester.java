package org.broadinstitute.hellbender.tools.variantdb.ingest.nextgen;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.*;
import org.broadinstitute.hellbender.tools.variantdb.ingest.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.ingest.IngestUtils;
import org.broadinstitute.hellbender.tools.variantdb.ingest.MetadataTsvCreator;
import org.broadinstitute.hellbender.utils.*;

import java.util.*;
import java.nio.file.Path;
import java.io.File;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Exome and Genome Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class VariantIngester extends VariantWalker {
    static final Logger logger = LogManager.getLogger(VariantIngester.class);

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

    @Argument(fullName = "output-path",
            shortName = "VPO",
            doc = "Path to the directory where the variants TSVs and positions expanded TSVs should be written")
    public GATKPathSpecifier parentOutputDirectory = null;
    public Path parentDirectory = null;

    @Argument(fullName = "ref-block-gq-to-ignore",
            shortName = "IG",
            doc = "Ref Block GQ band to ignore, bands of 10 e.g 0-9 get combined to 0, 20-29 get combined to 20",
            optional = true)
    public PetTsvCreator.GQStateEnum gqStateToIgnore = PetTsvCreator.GQStateEnum.SIXTY;

    @Argument(fullName = "sample-name-mapping",
            shortName = "SNM",
            doc = "Sample name to sample id mapping",
            optional = true)
    public File sampleMap;

    @Argument(fullName = "mode",
            shortName = "MO",
            doc = "Type of sample. Default is EXOMES. Valid options are EXOMES, GENOMES",
            optional = true)
    public CommonCode.ModeEnum mode = CommonCode.ModeEnum.EXOMES;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    private String refVersion = "37";


    @Override
    public boolean requiresIntervals() {
        return true; // TODO -- do I need to check the boolean flag on this?
    }

    @Override
    public void onTraversalStart() {

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);

        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        sampleName = IngestUtils.getSampleName(inputVCFHeader);
        sampleId = IngestUtils.getSampleId(sampleName, sampleMap);

        // Mod the sample directories
        int sampleDirectoryNumber = IngestUtils.getTableNumber(sampleId, IngestConstants.partitionPerTable);

        parentDirectory = parentOutputDirectory.toPath(); // TODO do we need this? More efficient way to do this?
        final Path sampleDirectoryPath = IngestUtils.createSampleDirectory(parentDirectory, sampleDirectoryNumber);
        metadataTsvCreator = new MetadataTsvCreator(sampleDirectoryPath);
        metadataTsvCreator.createRow(sampleName, sampleId, userIntervals, gqStateToIgnore);

        // To set up the missing positions
        SAMSequenceDictionary seqDictionary = getBestAvailableSequenceDictionary();
        userIntervals = intervalArgumentCollection.getIntervals(seqDictionary);

        final GenomeLocSortedSet genomeLocSortedSet = new GenomeLocSortedSet(new GenomeLocParser(seqDictionary));
        intervalArgumentGenomeLocSortedSet = GenomeLocSortedSet.createSetFromList(genomeLocSortedSet.getGenomeLocParser(), IntervalUtils.genomeLocsFromLocatables(genomeLocSortedSet.getGenomeLocParser(), intervalArgumentCollection.getIntervals(seqDictionary)));

        petTsvCreator = new PetTsvCreator(sampleName, sampleId, sampleDirectoryPath, seqDictionary, gqStateToIgnore);
        switch (mode) {
            case EXOMES:
                vetTsvCreator = new ExomeVetTsvCreator(sampleName, sampleId, sampleDirectoryPath);
                break;
            case GENOMES:
                vetTsvCreator = new GenomeVetTsvCreator(sampleName, sampleId, sampleDirectoryPath);
                break;
            case ARRAYS:
                throw new UserException.BadInput("To ingest Array data, use ArrayIngester tool.");
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
        if (!variant.isReferenceBlock() &&  variant.getAlternateAlleles().size() == 2 && variant.hasAlternateAllele(Allele.SPAN_DEL)){
            return;
        }


        // create VET output
        if (!variant.isReferenceBlock()) {
            vetTsvCreator.apply(variant, readsContext, referenceContext, featureContext);
        }
        petTsvCreator.apply(variant, intervalsToWrite);
    }

    @Override
    public Object onTraversalSuccess() {
        petTsvCreator.writeMissingIntervals(intervalArgumentGenomeLocSortedSet);
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
