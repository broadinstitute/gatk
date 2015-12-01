package org.broadinstitute.hellbender.tools.picard.vcf.filter;

import htsjdk.samtools.util.CollectionUtil;
import htsjdk.samtools.util.IOUtil;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderLineCount;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;

import java.io.File;
import java.util.*;

/**
 * Applies a set of hard filters to Variants and to Genotypes within a VCF.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Applies one or more hard filters to a VCF file to filter out genotypes and variants.",
        oneLineSummary = "Hard-filters a VCF file",
        programGroup = VariantProgramGroup.class
)
public final class FilterVcf extends PicardCommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, shortName= StandardArgumentDefinitions.INPUT_SHORT_NAME,
            doc="The input VCF file.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName=StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="The output VCF file.")
    public File OUTPUT;

    @Argument(doc="The minimum allele balance acceptable before filtering a site. Allele balance is calculated for heterozygotes as " +
            "the number of bases supporting the least-represented allele over the total number of base observations. Different heterozygote " +
            "genotypes at the same locus are measured independently. The locus is filtered if any allele balance is below the limit.")
    public double MIN_AB = 0.0d;

    @Argument(doc="The minimum sequencing depth supporting a genotype before the genotype will be filtered out.")
    public int MIN_DP = 0;

    @Argument(doc="The minimum genotype quality that must be achieved for a sample otherwise the genotype will be filtered out.")
    public int MIN_GQ = 0;

    @Argument(doc="The maximum phred scaled fisher strand value before a site will be filtered out.")
    public double MAX_FS = Double.MAX_VALUE;

    @Argument(doc="The minimum QD value to accept or otherwise filter out the variant.")
    public double MIN_QD = 0;

    /** Constructor to default to having index creation on. */
    public FilterVcf() { this.CREATE_INDEX = true; }

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsWritable(OUTPUT);

        final List<VariantFilter> variantFilters = CollectionUtil.makeList(new AlleleBalanceFilter(MIN_AB), new FisherStrandFilter(MAX_FS), new QdFilter(MIN_QD));
        final List<GenotypeFilter> genotypeFilters = CollectionUtil.makeList(new GenotypeQualityFilter(MIN_GQ), new DepthFilter(MIN_DP));
        try (final VCFFileReader in = new VCFFileReader(INPUT, false)) {
            final FilterApplyingVariantIterator iterator = new FilterApplyingVariantIterator(in.iterator(), variantFilters, genotypeFilters);

            try (final VariantContextWriter out = new VariantContextWriterBuilder().setOutputFile(OUTPUT).build()) {
                final VCFHeader header = in.getFileHeader();
                header.addMetaDataLine(new VCFFilterHeaderLine("AllGtsFiltered", "Site filtered out because all genotypes are filtered out."));
                header.addMetaDataLine(new VCFFormatHeaderLine("FT", VCFHeaderLineCount.UNBOUNDED, VCFHeaderLineType.String, "Genotype filters."));
                for (final VariantFilter filter : variantFilters) {
                    for (final VCFFilterHeaderLine line : filter.headerLines()) {
                        header.addMetaDataLine(line);
                    }
                }

                out.writeHeader(in.getFileHeader());

                while (iterator.hasNext()) {
                    out.add(iterator.next());
                }

            }
        }
        return null;
    }
}
