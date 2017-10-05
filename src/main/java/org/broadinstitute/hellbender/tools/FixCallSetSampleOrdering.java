package org.broadinstitute.hellbender.tools;

import com.google.common.collect.Lists;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLine;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.tools.genomicsdb.GenomicsDBImport;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.GenomicsDBTestUtils;
import picard.cmdline.programgroups.None;


import java.io.File;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

@CommandLineProgramProperties(summary = "fix the sample names in a vcf that ran through the GenotypeGVCFs pipeline when " +
        "the batch ordering bug was present, this will restore the correct sample names provided it is given the exact sample name mapping / vcf ordering" +
        "and batch sized that was used in the initial import",
        oneLineSummary = "fix sample names in a shuffled callset",
        programGroup = None.class,
        omitFromCommandLine = true)
public final class FixCallSetSampleOrdering extends VariantWalker {
    @Argument(fullName = GenomicsDBImport.SAMPLE_NAME_MAP_LONG_NAME,
            doc="the same sampleNameMap file which was used to generate the callset initially, it's important that this",
            optional = false)
    public String sampleNameMapPath;

    @Argument(fullName = GenomicsDBImport.BATCHSIZE_ARG_NAME,
            doc="the exact batch size that was used to generate the callset initially",
            minValue = 1,
            optional = false)
    public int batchSize;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="where to write the new unshuffled VCF",
            optional = false)
    public File output;

    private static VariantContextWriter writer;

    @Override
    public void onTraversalStart() {
        final List<String> sampleNamesOriginalOrdering = getSampleNamesInInputOrder(sampleNameMapPath);

        final VCFHeader originalHeader = getHeaderForVariants();
        final Set<VCFHeaderLine> originalHeaderLines = originalHeader.getMetaDataInInputOrder();

        final TreeSet<VCFHeaderLine> newHeaderLines = new TreeSet<>(originalHeaderLines);
        newHeaderLines.addAll(getDefaultToolVCFHeaderLines());
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

    private static List<String> getBatchSortedList(List<String> sampleNames, int batchSize){
        final List<List<String>> partition = Lists.partition(sampleNames, batchSize);
        return partition.stream()
                .flatMap( (List<String> batch) -> batch.stream().sorted() )
                .collect(Collectors.toList());
    }

}
