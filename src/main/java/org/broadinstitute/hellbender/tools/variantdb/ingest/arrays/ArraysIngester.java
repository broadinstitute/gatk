package org.broadinstitute.hellbender.tools.variantdb.ingest.arrays;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.variantdb.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.variantdb.ingest.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.ingest.IngestUtils;
import org.broadinstitute.hellbender.tools.variantdb.ingest.MetadataTsvCreator;

import java.io.File;
import java.nio.file.Path;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class ArraysIngester extends VariantWalker {
    static final Logger logger = LogManager.getLogger(ArraysIngester.class);

    private MetadataTsvCreator metadataTsvCreator;
    private RawArrayTsvCreator tsvCreator;

    private String sampleName;
    private String sampleId;

    // Inside the parent directory, a directory for each chromosome will be created, with a pet directory and vet directory in each one.
    // Each pet and vet directory will hold all of the pet and vet tsvs for each sample
    // A metadata directory will be created, with a metadata tsv for each sample

    @Argument(fullName = "output-path",
            shortName = "VPO",
            doc = "Path to the directory where the output TSVs should be written")
    public GATKPathSpecifier parentOutputDirectory = null;
    public Path parentDirectory = null;

    @Argument(fullName = "sample-name-mapping",
            shortName = "SNM",
            doc = "Sample name to sample id mapping",
            optional = true)
    public File sampleMap;

    @Argument(
            fullName = "ref-version",
            doc = "Remove this option!!!! only for ease of testing. Valid options are 37 or 38",
            optional = true)
    private String refVersion = "37";


    @Override
    public void onTraversalStart() {

        // Set reference version -- TODO remove this in the future, also, can we get ref version from the header?
        ChromosomeEnum.setRefVersion(refVersion);

        // Get sample name
        final VCFHeader inputVCFHeader = getHeaderForVariants();
        sampleName = IngestUtils.getSampleName(inputVCFHeader);
        sampleId = IngestUtils.getSampleId(sampleName, sampleMap);

        // Mod the sample directories
        int sampleDirectoryNumber = IngestUtils.getSampleDirectoryNumber(sampleId, IngestConstants.partitionPerTable);

        parentDirectory = parentOutputDirectory.toPath(); // TODO do we need this? More efficient way to do this?
        // If this sample set directory doesn't exist yet -- create it
        parentDirectory = parentOutputDirectory.toPath(); // TODO do we need this? More efficient way to do this?
        final Path sampleDirectoryPath = IngestUtils.createSampleDirectory(parentDirectory, sampleDirectoryNumber);
        metadataTsvCreator = new MetadataTsvCreator(sampleDirectoryPath);
        metadataTsvCreator.createRow(sampleName, sampleId);
        tsvCreator = new RawArrayTsvCreator(sampleName, sampleId, sampleDirectoryPath);
    }


    @Override
    public void apply(final VariantContext variant, final ReadsContext readsContext, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        tsvCreator.apply(variant, readsContext, referenceContext, featureContext);
    }

    @Override
    public Object onTraversalSuccess() {
        return 0;
    }

    @Override
    public void closeTool() {
        if (tsvCreator == null) {
            logger.info("tsvCreator is null when closing tool");
        }
        if (tsvCreator != null) {
            tsvCreator.closeTool();
        }
        if (metadataTsvCreator != null) {
            metadataTsvCreator.closeTool();
        }
    }
}
