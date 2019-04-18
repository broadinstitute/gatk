package org.broadinstitute.hellbender.tools.evoquer;

import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;

import java.util.List;

/**
 * Evoquer ("Evoker"):
 * Extract Variants Out of big QUERy.
 *
 * Connects to a BigQuery database and creates a VCF file from the given interval list based on the data in the database.
 *
 * The BigQuery database is assumed to have a schema and data that supports creation of {@link htsjdk.variant.variantcontext.VariantContext}
 * objects, and all data in this database are assumed to be based on the HG38 reference.
 *
 * An example invocation is:
 *
 *   ./gatk Evoquer --intervals chr2:10000-1741634 -R ~/references/Homo_sapiens_assembly38.fasta
 *
 * Created by jonn on 4/16/19.
 */
@CommandLineProgramProperties(
        summary = "(\"Evoker\") - Extract Variants Out of big QUERy.  Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        oneLineSummary = "Tool to extract VCF information from a BigQuery table containing genomic variant data (joing genotyping / GVCF-equivalent data).",
        programGroup = ShortVariantDiscoveryProgramGroup.class
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

    private EvoquerEngine evoquerEngine;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    /**
     * {@inheritDoc}
     *
     * {@link Evoquer} requires intervals.
     */
    @Override
    public boolean requiresIntervals() {
        return true;
    }

    /**
     * {@inheritDoc}
     *
     * {@link Evoquer} only requires a reference because the intervals need a sequence dictionary.
     * TODO: This should be fixed so that intervals can exist without a sequence dictionary.
     */
    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void onStartup() {
        super.onStartup();

        // Set up our EvoquerEngine:
        evoquerEngine = new EvoquerEngine();
    }

    @Override
    public void traverse() {
        final List<VariantContext> variants = evoquerEngine.evokeIntervals(getTraversalIntervals());

        logger.info( "Created the following variants:" );
        for ( final VariantContext variantContext : variants ) {
            logger.info( variantContext.toString() );
        }
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
