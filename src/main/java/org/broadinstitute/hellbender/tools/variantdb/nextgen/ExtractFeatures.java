package org.broadinstitute.hellbender.tools.variantdb.nextgen;

import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.tools.variantdb.CommonCode;
import org.broadinstitute.hellbender.tools.variantdb.SampleList;
import org.broadinstitute.hellbender.tools.variantdb.SchemaUtils;
import org.broadinstitute.hellbender.utils.bigquery.TableReference;

import java.util.HashSet;
import java.util.Set;

@CommandLineProgramProperties(
        summary = "(\"ExtractFeatures\") - Extract features data from BQ to train filtering model.",
        oneLineSummary = "Tool to extract variants out of big query to train filtering model",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public class ExtractFeatures extends ExtractTool {
    private static final Logger logger = LogManager.getLogger(ExtractFeatures.class);
    private ExtractFeaturesEngine engine;

    // TODO change this to vet_table_prefix!
    @Argument(
            fullName = "vet-table",
            doc = "Fully qualified name of the vet table",
            optional = false
    )
    protected String fqVetTable = null;

//    @Argument(
//            fullName = "pet-table",
//            doc = "Fully qualified name of the pet table where",
//            optional = false
//    )
//    protected String fqPetTable = null;

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

    @Override
    public boolean requiresIntervals() {
        return true; // TODO -- do I need to check the boolean flag on this?
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        TableReference sampleTableRef = new TableReference(sampleTableName, SchemaUtils.SAMPLE_FIELDS);

        VCFHeader header = CommonCode.generateVcfHeader(new HashSet<>(), reference.getSequenceDictionary());

        engine = new ExtractFeaturesEngine(
            projectID,
            vcfWriter,
            header,
            annotationEngine,
            reference,
            trainingSitesOnly,
            fqAltAlleleTable,
            fqVetTable,
            sampleTableRef,
            intervalArgumentCollection.getIntervals(getBestAvailableSequenceDictionary()),
            localSortMaxRecordsInRam,
            printDebugInformation,
            useBatchQueries,
            progressMeter);
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
