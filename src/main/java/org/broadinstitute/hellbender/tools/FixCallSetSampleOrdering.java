package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.GenomicsDBTestUtils;
import picard.cmdline.programgroups.None;


import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

@BetaFeature
@CommandLineProgramProperties(summary = "fix the sample names in a vcf that ran through the GenomicsDBImport when " +
        "the batch ordering bug was present, this will restore the correct sample names provided it is given the exact sample name mapping / vcf ordering" +
        "and batch sized that was used in the initial import",
        oneLineSummary = "fix sample names in a shuffled callset",
        programGroup = VariantProgramGroup.class,
        omitFromCommandLine = true)
public final class FixCallSetSampleOrdering extends VariantWalker {
    @Argument(fullName = GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME,
            doc="the same sampleNameMap file which was used to import the callset using GenomicsDBImport",
            optional = false)
    public String sampleNameMapPath;

    @Argument(fullName = GenomicsDBImport.BATCHSIZE_ARG_NAME,
            doc="the exact batch size that was used to import the callset using GenomicsDBImport",
            minValue = 0,
            optional = false)
    public int batchSize;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="where to write a reheadered version of the input VCF with the sample names in the correct order",
            optional = false)
    public File output;

    private VariantContextWriter writer;

    @Override
    public void onTraversalStart() {

        if (batchSize == 0) {
            throw new UserException("your callset is not affected by the bug if you ran with "+ GenomicsDBImport.BATCHSIZE_ARG_NAME +" 0");
        }

        final VCFHeader originalHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> originalHeaderLines = originalHeader.getMetaDataInInputOrder();
        final Set<VCFHeaderLine> newHeaderLines = new LinkedHashSet<>(originalHeaderLines);
        newHeaderLines.addAll(getDefaultToolVCFHeaderLines());

        final List<String> sampleNamesOriginalOrdering = getSampleNamesInInputOrder(sampleNameMapPath);
        if( sampleNamesOriginalOrdering.size() <= batchSize ){
            throw new UserException("you are not affected by the sample name ordering bug if your batch size is >= the number of samples in your callset. \n"
                                            + "batch size: " + batchSize + "\n"
                                            + "number of samples" + sampleNamesOriginalOrdering.size());
        }

        final List<String> batchSortedSampleNames = getBatchSortedList(sampleNamesOriginalOrdering, batchSize);
        final VCFHeader remappedHeader = new VCFHeader(newHeaderLines, batchSortedSampleNames);
        writer = createVCFWriter(output);
        writer.writeHeader(remappedHeader);
    }

    @Override
    public void apply(VariantContext variant, ReadsContext readsContext, ReferenceContext referenceContext, FeatureContext featureContext) {
        writer.add(variant);
    }

    @Override
    public void closeTool() {
        writer.close();
    }

    private static List<String> getSampleNamesInInputOrder(String sampleNameMapPath) {
        final Map<String, Path> stringPathMap = GenomicsDBImport.loadSampleNameMapFile(IOUtils.getPath(sampleNameMapPath));
        return new ArrayList<>(stringPathMap.keySet());
    }

    /**
     * Recreate the sort order that the buggy GenomicsDBImport used.  It takes the sample name in input order, partitions
     * them into batches, and then sorts within each batch.
     */
    private static List<String> getBatchSortedList(List<String> sampleNames, int batchSize){
        final List<List<String>> partition = Lists.partition(sampleNames, batchSize);
        return partition.stream()
                .flatMap( (List<String> batch) -> batch.stream().sorted() )
                .collect(Collectors.toList());
    }

}
