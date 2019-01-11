
package org.broadinstitute.hellbender.tools;

import htsjdk.samtools.util.IOUtil;
import htsjdk.tribble.*;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.util.TabixUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import picard.cmdline.programgroups.OtherProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.engine.ProgressMeter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.codecs.ProgressReportingDelegatingCodec;

import java.io.File;
import java.io.IOException;

/**
 * This tool creates an index file for the various kinds of feature-containing files supported by GATK (such as VCF
 * and BED files). An index allows querying features by a genomic interval.
 *
 * <h3>Usage example</h3>
 * <pre>
 * gatk IndexFeatureFile \
 *     -F cohort.vcf.gz
 * </pre>
 * This produces the corresponding index, cohort.vcf.gz.tbi.
 */

@CommandLineProgramProperties(
        summary = "Creates an index for a feature file, e.g. VCF or BED file.",
        oneLineSummary = "Creates an index for a feature file, e.g. VCF or BED file.",
        programGroup = OtherProgramGroup.class
)
@DocumentedFeature
public final class IndexFeatureFile extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger(IndexFeatureFile.class);

    @Argument(shortName = "F",
              fullName = "feature-file",
              doc = "Feature file (eg., VCF or BED file) to index. Must be in a tribble-supported format")
    public File featureFile;

    @Argument(shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
              fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
              doc = "The output index file. If missing, the tool will create an index file in the same directory " +
                     "as the input file.",
              optional = true)
    public File outputFile;

    public static final int OPTIMAL_GVCF_INDEX_BIN_SIZE = 128000;
    public static final String GVCF_FILE_EXTENSION = ".g.vcf";

    @Override
    protected Object doWork() {
        if (!featureFile.canRead()) {
            throw new UserException.CouldNotReadInputFile(featureFile);
        }

        // Get the right codec for the file to be indexed. This call will throw an appropriate exception
        // if featureFile is not in a supported format or is unreadable.
        final FeatureCodec<? extends Feature, ?> codec = new ProgressReportingDelegatingCodec<>(FeatureManager.getCodecForFile(featureFile), ProgressMeter.DEFAULT_SECONDS_BETWEEN_UPDATES);

        final Index index = createAppropriateIndexInMemory(codec);
        final File indexFile = determineFileName(index);

        try {
            index.write(indexFile);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile("Could not write index to file " + indexFile.getAbsolutePath(), e);
        }

        logger.info("Successfully wrote index to " + indexFile.getAbsolutePath());
        return indexFile.getAbsolutePath();
    }

    private File determineFileName(final Index index) {
        if (outputFile != null) {
            return outputFile;
        } else if (index instanceof TabixIndex) {
            return Tribble.tabixIndexFile(featureFile);
        } else {
            return Tribble.indexFile(featureFile);
        }
    }

    private Index createAppropriateIndexInMemory(final FeatureCodec<? extends Feature, ?> codec) {
        try {
            // For block-compression files, write a Tabix index
            if (IOUtil.hasBlockCompressedExtension(featureFile)) {
                // Creating tabix indices with a non standard extensions can cause problems so we disable it
                if (outputFile != null && !outputFile.getAbsolutePath().endsWith(TabixUtils.STANDARD_INDEX_EXTENSION)) {
                    throw new UserException("The index for " + featureFile + " must be written to a file with a \"" + TabixUtils.STANDARD_INDEX_EXTENSION + "\" extension");
                }

                // TODO: this could benefit from provided sequence dictionary from reference
                // TODO: this can be an optional parameter for the tool
                return IndexFactory.createIndex(featureFile, codec, IndexFactory.IndexType.TABIX, null);

            }
            // TODO: detection of GVCF files should not be file-extension-based. Need to come up with canonical
            // TODO: way of detecting GVCFs based on the contents (may require changes to the spec!)
            else if (featureFile.getName().endsWith(GVCF_FILE_EXTENSION)) {
                // Optimize GVCF indices for the use case of having a large number of GVCFs open simultaneously
                return IndexFactory.createLinearIndex(featureFile, codec, OPTIMAL_GVCF_INDEX_BIN_SIZE);
            } else {
                // Optimize indices for other kinds of files for seek time / querying
                return IndexFactory.createDynamicIndex(featureFile, codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
            }
        } catch (TribbleException e) {
            // Underlying cause here is usually a malformed file, but can also be things like
            // "codec does not support tabix"
            throw new UserException.CouldNotIndexFile(featureFile, e);
        }
    }
}
