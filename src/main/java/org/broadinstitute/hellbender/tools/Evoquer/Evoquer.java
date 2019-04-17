package org.broadinstitute.hellbender.tools.Evoquer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;

/**
 * Evoquer ("Evoker"):
 * Extract Variants Out of big QUERy.
 *
 * Connects to a BigQuery database and creates a VCF file from the given interval list based on the data in the database.
 *
 * The BigQuery database is assumed to have a schema and data that supports creation of {@link htsjdk.variant.variantcontext.VariantContext}
 * objects, and all data in this database are assumed to be based on the HG38 reference.
 *
 * Created by jonn on 4/16/19.
 */
@CommandLineProgramProperties(
        summary = "(\"Evoker\") - Extract Variants Out of big QUERy.  Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        oneLineSummary = "Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        programGroup = ShortVariantDiscoveryProgramGroup.class,
        omitFromCommandLine = false
)
@DocumentedFeature
public class Evoquer extends GATKTool {
    private static final Logger logger = LogManager.getLogger(Evoquer.class);

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
