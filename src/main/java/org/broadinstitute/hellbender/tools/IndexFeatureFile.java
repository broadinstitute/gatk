package org.broadinstitute.hellbender.tools;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.Tribble;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.index.tabix.TabixIndex;
import htsjdk.tribble.util.LittleEndianOutputStream;
import htsjdk.variant.vcf.VCF3Codec;
import htsjdk.variant.vcf.VCFCodec;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureManager;
import org.broadinstitute.hellbender.exceptions.UserException;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Tool to create an appropriate index file for the various kinds of Feature-containing files
 * we support. These include VCF, BCF, and BED files.
 *
 * Such files must have an index in order to be queried by interval.
 */
@CommandLineProgramProperties(
        usage = "Creates indices for Feature-containing files, such as VCF, BCF, and BED files",
        usageShort = "Creates indices for Feature-containing files",
        programGroup = VariantProgramGroup.class
)
public final class IndexFeatureFile extends CommandLineProgram {
    private static final Logger logger = LogManager.getLogger(IndexFeatureFile.class);

    @Argument(shortName = "F", fullName = "feature_file", doc = "Feature file (eg., VCF/BCF/etc. file) to index. Must be in a tribble-supported format")
    private File featureFile;

    public static final int OPTIMAL_GVCF_INDEX_BIN_SIZE = 128000;
    public static final String GVCF_FILE_EXTENSION = ".gvcf";
    public static final String TABIX_INDEX_EXTENSION = ".tbi";

    @Override
    protected Object doWork() {
        if ( ! featureFile.canRead() ) {
            throw new UserException.CouldNotReadInputFile(featureFile);
        }

        // Get the right codec for the file to be indexed. This call will throw an appropriate exception
        // if featureFile is not in a supported format or is unreadable.
        FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(featureFile);

        // Create appropriate index type in memory, and determine the name of the file we'll write the index to
        final Index index = createAppropriateIndexInMemory(featureFile, codec);
        final File indexFile = index instanceof TabixIndex ? new File(featureFile.getAbsolutePath() + TABIX_INDEX_EXTENSION) :
                                                             Tribble.indexFile(featureFile);

        // Write the index to disk
        try ( LittleEndianOutputStream stream = new LittleEndianOutputStream(new FileOutputStream(indexFile)) ) {
            index.write(stream);
        }
        catch ( IOException e ) {
            throw new UserException.CouldNotCreateOutputFile("Could not write index to file " + indexFile.getAbsolutePath(), e);
        }

        logger.info("Successfully wrote index to " + indexFile.getAbsolutePath());
        return 0;
    }

    private Index createAppropriateIndexInMemory( final File featureFile, final FeatureCodec<? extends Feature, ?> codec ) {
        // For block-compression VCF files, write a Tabix index
        if ( AbstractFeatureReader.hasBlockCompressedExtension(featureFile) ) {
            if ( isVCFCodec(codec) ) {
                return IndexFactory.createTabixIndex(featureFile, codec, TabixFormat.VCF, null);
            }
            else {
                throw new UserException("This tool only supports indexing of block-compressed files when they are in VCF format");
            }
        }
        // TODO: detection of GVCF files should not be file-extension-based. Need to come up with canonical
        // TODO: way of detecting GVCFs based on the contents (may require changes to the spec!)
        else if ( featureFile.getName().endsWith(GVCF_FILE_EXTENSION) ) {
            // Optimize GVCF indices for the use case of having a large number of GVCFs open simultaneously
            return IndexFactory.createLinearIndex(featureFile, codec, OPTIMAL_GVCF_INDEX_BIN_SIZE);
        }
        else {
            // Optimize indices for other kinds of files for seek time / querying
            return IndexFactory.createDynamicIndex(featureFile, codec, IndexFactory.IndexBalanceApproach.FOR_SEEK_TIME);
        }
    }

    private boolean isVCFCodec( final FeatureCodec<? extends Feature, ?> codec ) {
        return codec.getClass() == VCFCodec.class || codec.getClass() == VCF3Codec.class;
    }
}
