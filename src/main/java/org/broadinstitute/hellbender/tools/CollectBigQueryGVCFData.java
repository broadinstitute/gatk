package org.broadinstitute.hellbender.tools;

import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;

/**
 * Connects to a BigQuery database and creates a VCF file from the given interval list based on the data in the database.
 *
 * The BigQuery database is assumed to have a schema and data that supports creation of {@link htsjdk.variant.variantcontext.VariantContext}
 * objects.
 *
 * Created by jonn on 4/16/19.
 */
@CommandLineProgramProperties(
        summary = "Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        oneLineSummary = "Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = false
)
@DocumentedFeature
public class CollectBigQueryGVCFData extends GATKTool {

    //==================================================================================================================
    // Public Static Members:

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public void traverse() {

    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
