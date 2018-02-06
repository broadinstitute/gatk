package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.TabixUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ShortVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by sam on 2/5/18.
 */
@CommandLineProgramProperties(
        summary = "Apply tranche filtering",
        oneLineSummary = "Apply tranche filtering",
        programGroup = ShortVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public class NeuralNetTranches extends CommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc = "Input VCF file")
    private String inputVcf = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private String outputVcf = null;

    @Argument(fullName = "snpTruthVcf", shortName = "stv", doc = "Input file of known common SNP sites.")
    private String snpTruthVcf = null;

    @Argument(fullName = "indelTruthVcf", shortName = "itv", doc = "Input file of known common INDEL sites.")
    private String indelTruthVcf = null;

    @Argument(fullName = "scoreKey", shortName = "sk", doc = "Score key.")
    private String scoreKey = "CNN_1D";

    @Argument(fullName = "samples", shortName = "s", doc = "Maximum number of truth VCF sites to check.")
    private int samples = 1200;


    @Argument(fullName="tranche",
            shortName="t",
            doc="The levels of truth sensitivity at which to slice the data. (in percents, that is 1.0 for 1 percent)",
            optional=true)
    private List<Double> tranches = new ArrayList<Double>(Arrays.asList(99.9, 99.0, 90.0));

    // Start the Python executor. This does not actually start the Python process, but fails if python can't be located
    final PythonScriptExecutor pythonExecutor = new PythonScriptExecutor(true);
    private File tempFile, tempFileIdx;
    private String trancheString;

    @Override
    protected void onStartup() {
        /* check for successful import of libraries */
        PythonScriptExecutor.checkPythonEnvironmentForPackage("vqsr_cnn");
        trancheString = " --tranches "+tranches.toString().substring(1, tranches.toString().length()-1).replace(",", "");

        try {
            final String idxExt = ".tbi";
            tempFile = File.createTempFile(outputVcf, "_temp.vcf.gz");
            tempFileIdx = new File(tempFile.getAbsolutePath()+idxExt);
        } catch (IOException e) {
            throw new GATKException("Error when creating temp files.", e);
        }
    }

    @Override
    protected Object doWork() {
        final Resource pythonScriptResource = new Resource("tranches.py", org.broadinstitute.hellbender.tools.walkers.vqsr.NeuralNetTranches.class);
        final List<String> snpArguments = new ArrayList<>(Arrays.asList(
                "--mode", "write_snp_tranches",
                "--input_vcf", inputVcf,
                "--train_vcf", snpTruthVcf,
                "--score_keys", scoreKey,
                trancheString,
                "--samples", Integer.toString(samples),
                "--output_vcf", tempFile.getAbsolutePath()));

        logger.info("SNP Args are:"+ Arrays.toString(snpArguments.toArray()));
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                snpArguments
        );

        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(tempFile);
        createAppropriateIndexInMemory(codec, tempFile, tempFileIdx);
        final Index index = createAppropriateIndexInMemory(codec, tempFile, tempFileIdx);

        try {
            index.write(tempFileIdx);
        } catch (final IOException e) {
            throw new GATKException(String.format("Could not write temporary index to file: %s", tempFileIdx.getAbsolutePath()), e);
        }

        final List<String> indelArguments = new ArrayList<>(Arrays.asList(
                "--mode", "write_indel_tranches",
                "--input_vcf", tempFile.getAbsolutePath(),
                "--train_vcf", indelTruthVcf,
                "--score_keys", scoreKey,
                trancheString,
                "--samples", Integer.toString(samples),
                "--output_vcf", outputVcf));

        logger.info("Did SNP filtering, now INDELs. Arguments are:"+ Arrays.toString(indelArguments.toArray()));
        final boolean pythonReturnCode2 = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                indelArguments
        );
        return pythonReturnCode && pythonReturnCode2;
    }

    // Stolen and adapted from IndexFeatureFile.java
    private Index createAppropriateIndexInMemory(final FeatureCodec<? extends Feature, ?> codec, File featureFile, File indexFile) {
        try {
            // For block-compression files, write a Tabix index
            if (AbstractFeatureReader.hasBlockCompressedExtension(featureFile)) {
                // Creating tabix indices with a non standard extensions can cause problems so we disable it
                if (!indexFile.getAbsolutePath().endsWith(TabixUtils.STANDARD_INDEX_EXTENSION)) {
                    throw new UserException("The index for " + featureFile + " must be written to a file with a \"" + TabixUtils.STANDARD_INDEX_EXTENSION + "\" extension");
                }

                return IndexFactory.createIndex(featureFile, codec, IndexFactory.IndexType.TABIX, null);

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
