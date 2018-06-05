package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.Feature;
import htsjdk.tribble.FeatureCodec;
import htsjdk.tribble.TribbleException;
import htsjdk.tribble.index.Index;
import htsjdk.tribble.index.IndexFactory;
import htsjdk.tribble.util.TabixUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import picard.cmdline.programgroups.VariantFilteringProgramGroup;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Apply tranche filtering to VCF based on scores from the INFO field.
 *
 * <h3>Inputs</h3>
 * <ul>
 *      <li>The input variants to tranche filter.</li>
 *      <li>snp-truth-vcf A VCF containing known common SNP sites</li>
 *      <li>indel-truth-vcf A VCF containing known and common INDEL sites.</li>
 *      <li>info-key The key from the INFO field of the VCF which contains the values that will be used to filter.</li>
 *      <li>tranche List of percent sensitivities to the known sites at which we will filter.  Must be between 0 and 100.</li>
 * </ul>
 *
 * <h3>Outputs</h3>
 * <ul>
 * <li>A tranche filtered VCF.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 *
 * <h4>Apply tranche filters based on CNN_1D scores</h4>
 * <pre>
 * gatk FilterVariantTranches \
 *   -V input.vcf.gz \
 *   --snp-truth-vcf hapmap.vcf \
 *   --indel-truth-vcf mills.vcf \
 *   --info-key CNN_1D \
 *   --tranche 99.9 --tranche 99.0 --tranche 95 \
 *   --max-sites 8000 \
 *   -O filtered.vcf
 * </pre>
 *
 */
@CommandLineProgramProperties(
        summary = "Apply tranche filtering based on a truth VCF of known common sites of variation and score from VCF INFO field",
        oneLineSummary = "Apply tranche filtering",
        programGroup = VariantFilteringProgramGroup.class
)
@DocumentedFeature
@ExperimentalFeature
public class FilterVariantTranches extends CommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME,
            doc = "Input VCF file")
    private String inputVcf = null;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc = "Output VCF file")
    private String outputVcf = null;

    @Argument(fullName = "snp-truth-vcf", shortName = "snp-truth-vcf", doc = "Input file of known common SNP sites.")
    private String snpTruthVcf = null;

    @Argument(fullName = "indel-truth-vcf", shortName = "indel-truth-vcf", doc = "Input file of known common INDEL sites.")
    private String indelTruthVcf = null;

    @Argument(fullName = "info-key", shortName = "info-key", doc = "The key must be in the INFO field of the input VCF.")
    private String infoKey = GATKVCFConstants.CNN_1D_KEY;

    @Argument(fullName = "max-sites", shortName = "max-sites", doc = "Maximum number of truth VCF sites to check.")
    private int maxSites = 1200;

    @Argument(fullName="tranche",
            shortName="t",
            doc="The levels of truth sensitivity at which to slice the data. (in percents, i.e. 99.9 for 99.9 percent and 1.0 for 1 percent)",
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
        trancheString = " --tranches " + tranches.stream().map(d -> Double.toString(d)).collect(Collectors.joining(" "));
        try {
            final String idxExt = ".tbi";
            tempFile = File.createTempFile(outputVcf, "_temp.vcf.gz");
            tempFile.deleteOnExit();
            tempFileIdx = new File(tempFile.getAbsolutePath()+idxExt);
            tempFileIdx.deleteOnExit();
        } catch (IOException e) {
            throw new GATKException("Error when creating temp files.", e);
        }
    }

    @Override
    protected Object doWork() {
        final Resource pythonScriptResource = new Resource("tranches.py", FilterVariantTranches.class);
        final List<String> snpArguments = new ArrayList<>(Arrays.asList(
                "--mode", "write_snp_tranches",
                "--input_vcf", inputVcf,
                "--train_vcf", snpTruthVcf,
                "--score_keys", infoKey,
                trancheString,
                "--samples", Integer.toString(maxSites),
                "--output_vcf", tempFile.getAbsolutePath()));

        logger.info("SNP Args are:"+ Arrays.toString(snpArguments.toArray()));
        final boolean pythonReturnCode = pythonExecutor.executeScript(
                pythonScriptResource,
                null,
                snpArguments
        );

        final FeatureCodec<? extends Feature, ?> codec = FeatureManager.getCodecForFile(tempFile);
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
                "--score_keys", infoKey,
                trancheString,
                "--samples", Integer.toString(maxSites),
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
