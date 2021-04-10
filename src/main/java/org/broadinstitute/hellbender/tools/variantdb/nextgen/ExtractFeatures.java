package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SampleList;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.util.HashSet;
import java.util.List;

@CommandLineProgramProperties(
        summary = "(\"ExtractFeatures\") - Extract features data from BQ to train filtering model.",
        oneLineSummary = "Tool to extract variants out of big query to train filtering model",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractFeatures extends ExtractTool {
    private static final Logger logger = LogManager.getLogger(ExtractFeatures.class);
    private ExtractFeaturesEngine engine;

    @Argument(
            fullName = "alt-allele-table",
            doc = "Fully qualified name of the table where the alternate allele info is",
            optional = false
    )
    protected String fqAltAlleleTable = null;

    @Argument(
            fullName = "training-sites-only",
            doc = "Whether the extract is for training the model (a limited number of sites) or applying the model. Default is applying the model.",
            optional = true
    )
    protected boolean trainingSitesOnly = false;

    @Argument(
        fullName = "use-batch-queries",
        doc = "If true, use batch (rather than interactive) priority queries in BigQuery",
        optional = true)
    protected boolean useBatchQueries = true;

    @Argument(
        fullName = "hq-genotype-gq-threshold",
        doc = "GQ threshold defining a high quality genotype",
        optional = true)
    protected int hqGenotypeGQThreshold = 20;

    @Argument(
        fullName = "hq-genotype-depth-threshold",
        doc = "Depth threshold defining a high quality genotype",
        optional = true)
    protected int hqGenotypeDepthThreshold = 10;

    @Argument(
        fullName = "hq-genotype-ab-threshold",
        doc = "Ab threshold defining a high quality genotype",
        optional = true)
    protected double hqGenotypeABThreshold = 0.2;

    @Override
    public boolean requiresIntervals() {
        return false;
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        TableReference sampleTableRef = new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS);
        SampleList sampleList = new SampleList(sampleTableName, sampleFileName, projectID, printDebugInformation);

        VCFHeader header = CommonCode.generateVcfHeader(new HashSet<>(), reference.getSequenceDictionary());

        final List<SimpleInterval> traversalIntervals = getTraversalIntervals();

        if (minLocation == null && maxLocation == null && hasUserSuppliedIntervals()) {
            final SimpleInterval firstInterval = traversalIntervals.get(0);
            final SimpleInterval lastInterval = traversalIntervals.get(traversalIntervals.size() - 1);

            minLocation = SchemaUtils.encodeLocation(firstInterval.getContig(), firstInterval.getStart());
            maxLocation = SchemaUtils.encodeLocation(lastInterval.getContig(), lastInterval.getEnd());
        } else if ((minLocation != null || maxLocation != null) && hasUserSuppliedIntervals()) {
            throw new UserException("min-location and max-location should not be used together with intervals (-L).");
        }

        engine = new ExtractFeaturesEngine(
            projectID,
            vcfWriter,
            header,
            annotationEngine,
            reference,
            trainingSitesOnly,
            fqAltAlleleTable,
            sampleTableRef,
            traversalIntervals,
            minLocation,
            maxLocation,
            localSortMaxRecordsInRam,
            printDebugInformation,
            useBatchQueries,
            progressMeter,
            sampleList.size(),
            hqGenotypeGQThreshold,
            hqGenotypeDepthThreshold,
            hqGenotypeABThreshold);

        vcfWriter.writeHeader(header);
}

    @Override
    // maybe think about creating a BigQuery Row walker?
    public void traverse() {
        progressMeter.setRecordsBetweenTimeChecks(100L);

        engine.traverse();
    }

    @Override
    protected void onShutdown() {
        super.onShutdown();

        // Close up our writer if we have to:
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }
    }

}
