package org.broadinstitute.hellbender.tools.variantdb.arrays;

import htsjdk.samtools.util.RuntimeIOException;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.variantdb.ChromosomeEnum;
import org.broadinstitute.hellbender.tools.variantdb.IngestConstants;
import org.broadinstitute.hellbender.tools.variantdb.IngestUtils;
import org.broadinstitute.hellbender.tools.variantdb.arrays.tables.ProbeInfo;

import java.io.File;
import java.util.Map;

/**
 * Ingest variant walker
 */
@CommandLineProgramProperties(
        summary = "Ingest tool for the Joint Genotyping in Big Query project",
        oneLineSummary = "Ingest tool for BQJG",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = true
)
public final class CreateArrayIngestFiles extends VariantWalker {
    static final Logger logger = LogManager.getLogger(CreateArrayIngestFiles.class);

    private ArraySampleTsvCreator sampleTsvCreator;
    private RawArrayTsvCreator tsvCreator;

    private String sampleName;
    private String sampleId;

    @Argument(
            fullName = "metrics-file",
            shortName = "QCF",
            doc = "Filepath to picard metrics file",
            optional = true)  // TODO change this to false for release
    private String metricsFilePath = null;

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

    @Argument(fullName = "probe-info-table",
            shortName = "PIT",
            doc = "Fully qualified table name for the probe info table",
            optional = true)
    public String probeFQTablename;

    @Argument(
        fullName = "probe-info-file",
        shortName = "PIF",
        doc = "Filepath to CSV export of probe-info table",
        optional = true)
    private String probeCsvFile = null;


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
    public void onTraversalStart() {
        //set up output directory
        if (!outputDir.exists() && !outputDir.mkdir()) {
            throw new RuntimeIOException("Unable to create directory: " + outputDir.getAbsolutePath());
        }

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

        sampleTsvCreator = new ArraySampleTsvCreator(metricsFilePath);
        sampleTsvCreator.createRow(sampleName, sampleId, tableNumberPrefix, outputDir);

        Map<String, ProbeInfo> probeNameMap;
        if (probeCsvFile == null) {
            probeNameMap = ExtractCohortBQ.getProbeNameMap(probeFQTablename, false);
        } else {
            probeNameMap = ProbeInfo.getProbeNameMap(probeCsvFile);
        }

        // Set reference version
        ChromosomeEnum.setRefVersion(refVersion);

        tsvCreator = new RawArrayTsvCreator(sampleName, sampleId, tableNumberPrefix, probeNameMap, outputDir);
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
        if (sampleTsvCreator != null) {
            sampleTsvCreator.closeTool();
        }
    }
}
