package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.MergingIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.VariantContextComparator;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;
import java.util.stream.StreamSupport;

/**
 * Mutect2 and its filtering engine optionally use a panel of normals to reduce false positive calls.  The idea is to
 * run Mutect2 -- without filtering -- on a large collection of normal samples.  Any resulting call is most likely a
 * sequencing artifact or a germline variant, neither of which are true somatic calls.  This tool takes these separate
 * normal vcfs and creates a single vcf of false positive variants.
 *
 * Use:
 *
 * 1) Run Mutect without filtering on a large number of normal samples (note that the argument "-tumor normal1SampleName"
 * is NOT a typo -- the idea is to treat a normal as if it were a tumor in tumor-only mode)
 *  java -jar gatk.jar Mutect2 -R reference.fasta -I normal1.bam -tumor normal1SampleName -o output1.vcf
 *  java -jar gatk.jar Mutect2 -R reference.fasta -I normal2.bam -tumor normal1SampleName -o output2.vcf
 *  . . .
 *
 *  2) Create one or more files ending in ".list" with all the paths to the vcfs in step 1.  The following are the
 *     contents of vcfs.list:
 *  output1.vcf
 *  output2.vcf
 *  . . .
 *
 *  3) Run this tool
 *  java jar gatk.jar CreateSomaticPanelOfNormals -vcfs vcfs.list -O pon.vcf
 *
 *  4) Use the panel of normals vcf to improve Mutect calls
 *  java -jar gatk.jar -R reference.fasta -I tumor.bam -tumor tumorSampleName \
 *    [-I matchedNormal.bam -normal matchedNormalSampleName] \
 *    -PON pon.vcf
 *
 *  Note that multiple .list files may be passed in by specifying the -vcfs option multiple times. It's also possible
 *  to just provide the VCFs explicitly on the command line, rather than via .list files.
 *
 * Created by David Benjamin on 2/17/17.
 */
@CommandLineProgramProperties(
        summary = "Make a somatic panel of normals",
        oneLineSummary = "Make a somatic panel of normals",
        programGroup = VariantProgramGroup.class
)
public class CreateSomaticPanelOfNormals extends CommandLineProgram {

    public static final String INPUT_VCFS_LIST_LONG_NAME = "vcfsListFile";
    public static final String INPUT_VCFS_LIST_SHORT_NAME = "vcfs";

    /**
     * The VCFs can be input as either one or more .list file(s) containing one VCF per line, or VCFs can be
     * specified explicitly on the command line.
     */
    @Argument(fullName = INPUT_VCFS_LIST_LONG_NAME,
            shortName = INPUT_VCFS_LIST_SHORT_NAME,
            doc="VCFs for samples to include. May be specified either one at a time, or as one or more .list file containing multiple VCFs, one per line.", optional = false)
    private Set<File> vcfs = new LinkedHashSet<>(0);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="Output vcf", optional = false)
    private File outputVcf = null;

    public Object doWork() {
        final List<File> inputVcfs = new ArrayList<>(vcfs);
        final Collection<CloseableIterator<VariantContext>> iterators = new ArrayList<>(inputVcfs.size());
        final Collection<VCFHeader> headers = new HashSet<>(inputVcfs.size());
        final VCFHeader headerOfFirstVcf = new VCFFileReader(inputVcfs.get(0), false).getFileHeader();
        final SAMSequenceDictionary sequenceDictionary = headerOfFirstVcf.getSequenceDictionary();
        final VariantContextComparator comparator = headerOfFirstVcf.getVCFRecordComparator();


        for (final File vcf : inputVcfs) {
            final VCFFileReader reader = new VCFFileReader(vcf, false);
            iterators.add(reader.iterator());
            final VCFHeader header = reader.getFileHeader();
            Utils.validateArg(comparator.isCompatible(header.getContigLines()), () -> vcf.getAbsolutePath() + " has incompatible contigs.");
            headers.add(header);
        }

        final VariantContextWriter writer = GATKVariantContextUtils.createVCFWriter(outputVcf, sequenceDictionary, false, Options.INDEX_ON_THE_FLY);
        writer.writeHeader(new VCFHeader(VCFUtils.smartMergeHeaders(headers, false)));

        final MergingIterator<VariantContext> mergingIterator = new MergingIterator<>(comparator, iterators);
        SimpleInterval currentPosition = new SimpleInterval("FAKE", 1, 1);
        final List<VariantContext> variantsAtThisPosition = new ArrayList<>(20);
        while (mergingIterator.hasNext()) {
            final VariantContext vc = mergingIterator.next();
            if (!currentPosition.overlaps(vc)) {
                processVariantsAtSamePosition(variantsAtThisPosition, writer);
                variantsAtThisPosition.clear();
                currentPosition = new SimpleInterval(vc.getContig(), vc.getStart(), vc.getStart());
            }
            variantsAtThisPosition.add(vc);
        }
        mergingIterator.close();
        writer.close();

        return "SUCCESS";
    }

    //TODO: this is the old Mutect behavior that just looks for multiple hits
    //TODO: we should refine this
    private static void processVariantsAtSamePosition(final List<VariantContext> variants, final VariantContextWriter writer) {
        if (variants.size() > 1){
            final VariantContext mergedVc = AssemblyBasedCallerUtils.makeMergedVariantContext(variants);
            final VariantContext outputVc = new VariantContextBuilder()
                    .source(mergedVc.getSource())
                    .loc(mergedVc.getContig(), mergedVc.getStart(), mergedVc.getEnd())
                    .alleles(mergedVc.getAlleles())
                    .make();
            writer.add(outputVc);
        }
    }
}
