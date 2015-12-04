package org.broadinstitute.hellbender.tools.picard.vcf;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.liftover.LiftOver;
import htsjdk.samtools.reference.ReferenceSequenceFileWalker;
import htsjdk.samtools.util.*;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFFileReader;
import htsjdk.variant.vcf.VCFFilterHeaderLine;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFRecordCodec;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.PicardCommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.utils.runtime.ProgressLogger;

import java.io.File;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Tool for lifting over a VCF to another genome build and producing a properly header'd,
 * sorted and indexed VCF in one go.
 *
 * @author Tim Fennell
 */
@CommandLineProgramProperties(
        summary = "Lifts a VCF over from one genome build to another using UCSC liftover. The output file will be sorted " +
                "and indexed. Records may be rejected because they cannot be lifted over or because post-liftover the " +
                "reference allele mismatches the target genome build.  Rejected records will be emitted with filters " +
                "to the REJECT file, on the source genome.",
        oneLineSummary = "Lifts a VCF between genome builds",
        programGroup = VariantProgramGroup.class
)
public final class LiftOverVcf extends PicardCommandLineProgram {
    @Argument(fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME, 
            doc = "The input VCF file to be lifted over.")
    public File INPUT;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, 
            doc = "The output location to write the lifted over VCF to.")
    public File OUTPUT;

    @Argument(shortName = "C", doc = "The liftover chain file. See https://genome.ucsc.edu/goldenPath/help/chain.html for a description" +
            " of chain files.  See http://hgdownload.soe.ucsc.edu/downloads.html#terms for where to download chain files.")
    public File CHAIN;

    @Argument(doc = "File to which to write rejected records.")
    public File REJECT;

    /** Filter name to use when a target cannot be lifted over. */
    public static final String FILTER_CANNOT_LIFTOVER = "FailedLiftover";

    /** Filter name to use when a target is lifted over, but the reference allele doens't match the new reference. */
    public static final String FILTER_MISMATCHING_REF_ALLELE = "MismatchedRefAllele";

    /** Filters to be added to the REJECT file. */
    private static final List<VCFFilterHeaderLine> FILTERS = CollectionUtil.makeList(
            new VCFFilterHeaderLine(FILTER_CANNOT_LIFTOVER, "Variant could not be lifted between genome builds."),
            new VCFFilterHeaderLine(FILTER_MISMATCHING_REF_ALLELE, "Reference allele does not match reference genome sequence after liftover.")
    );

    @Override
    protected Object doWork() {
        IOUtil.assertFileIsReadable(INPUT);
        IOUtil.assertFileIsReadable(REFERENCE_SEQUENCE);
        IOUtil.assertFileIsReadable(CHAIN);
        IOUtil.assertFileIsWritable(OUTPUT);
        IOUtil.assertFileIsWritable(REJECT);

        ////////////////////////////////////////////////////////////////////////
        // Setup the inputs
        ////////////////////////////////////////////////////////////////////////
        final LiftOver liftOver = new LiftOver(CHAIN);
        final VCFFileReader in = new VCFFileReader(INPUT, false);

        logger.info("Loading up the target reference genome.");
        final ReferenceSequenceFileWalker walker = new ReferenceSequenceFileWalker(REFERENCE_SEQUENCE);
        final Map<String,byte[]> refSeqs = new HashMap<>();
        for (final SAMSequenceRecord rec: walker.getSequenceDictionary().getSequences()) {
            refSeqs.put(rec.getSequenceName(), walker.get(rec.getSequenceIndex()).getBases());
        }
        CloserUtil.close(walker);


        ////////////////////////////////////////////////////////////////////////
        // Setup the outputs
        ////////////////////////////////////////////////////////////////////////
        final VCFHeader inHeader = in.getFileHeader();
        final VCFHeader outHeader = new VCFHeader(inHeader);
        outHeader.setSequenceDictionary(walker.getSequenceDictionary());
        final VariantContextWriter out = new VariantContextWriterBuilder().setOption(Options.INDEX_ON_THE_FLY)
                .setOutputFile(OUTPUT).setReferenceDictionary(walker.getSequenceDictionary()).build();
        out.writeHeader(outHeader);

        final VariantContextWriter rejects = new VariantContextWriterBuilder().setOutputFile(REJECT).unsetOption(Options.INDEX_ON_THE_FLY).build();
        final VCFHeader rejectHeader = new VCFHeader(in.getFileHeader());
        for (final VCFFilterHeaderLine line : FILTERS) rejectHeader.addMetaDataLine(line);
        rejects.writeHeader(rejectHeader);


        ////////////////////////////////////////////////////////////////////////
        // Read the input VCF, lift the records over and write to the sorting
        // collection.
        ////////////////////////////////////////////////////////////////////////
        long failedLiftover = 0, failedAlleleCheck = 0, total = 0;
        logger.info("Lifting variants over and sorting.");

        final SortingCollection<VariantContext> sorter = SortingCollection.newInstance(VariantContext.class,
                new VCFRecordCodec(outHeader),
                outHeader.getVCFRecordComparator(),
                MAX_RECORDS_IN_RAM,
                TMP_DIR);

        ProgressLogger progress = new ProgressLogger(logger, 1000000, "read");

        for (final VariantContext ctx : in) {
            ++total;
            final Interval source = new Interval(ctx.getContig(), ctx.getStart(), ctx.getEnd(), false, ctx.getContig() + ":" + ctx.getStart() + "-" + ctx.getEnd());
            final Interval target = liftOver.liftOver(source, 1.0);

            if (target == null) {
                rejects.add(new VariantContextBuilder(ctx).filter(FILTER_CANNOT_LIFTOVER).make());
                failedLiftover++;
            }
            else {
                // Fix the alleles if we went from positive to negative strand
                final List<Allele> alleles = new ArrayList<>();
                for (final Allele oldAllele : ctx.getAlleles()) {
                    if (target.isPositiveStrand() || oldAllele.isSymbolic()) {
                        alleles.add(oldAllele);
                    }
                    else {
                        alleles.add(Allele.create(SequenceUtil.reverseComplement(oldAllele.getBaseString()), oldAllele.isReference()));
                    }
                }

                // Build the new variant context
                final VariantContextBuilder builder = new VariantContextBuilder(
                        ctx.getSource(),
                        target.getContig(),
                        target.getStart(),
                        target.getEnd(),
                        alleles);

                builder.id(ctx.getID());
                builder.attributes(ctx.getAttributes());
                builder.genotypes(ctx.getGenotypes());
                builder.filters(ctx.getFilters());
                builder.log10PError(ctx.getLog10PError());

                // Check that the reference allele still agrees with the reference sequence
                boolean mismatchesReference = false;
                for (final Allele allele : builder.getAlleles()) {
                    if (allele.isReference()) {
                        final byte[] ref = refSeqs.get(target.getContig());
                        final String refString = StringUtil.bytesToString(ref, target.getStart() - 1, target.length());

                        if (!refString.equalsIgnoreCase(allele.getBaseString())) {
                            mismatchesReference = true;
                        }

                        break;
                    }
                }

                if (mismatchesReference) {
                    rejects.add(new VariantContextBuilder(ctx).filter(FILTER_MISMATCHING_REF_ALLELE).make());
                    failedAlleleCheck++;
                }
                else {
                    sorter.add(builder.make());
                }
            }

            progress.record(ctx.getContig(), ctx.getStart());
        }

        final NumberFormat pfmt = new DecimalFormat("0.0000%");
        final String pct = pfmt.format((failedLiftover + failedAlleleCheck) / (double) total);
        logger.info("Processed ", total, " variants.");
        logger.info(Long.toString(failedLiftover), " variants failed to liftover.");
        logger.info(Long.toString(failedAlleleCheck), " variants lifted over but had mismatching reference alleles after lift over.");
        logger.info(pct, " of variants were not successfully lifted over and written to the output.");

        rejects.close();
        in.close();

        ////////////////////////////////////////////////////////////////////////
        // Write the sorted outputs to the final output file
        ////////////////////////////////////////////////////////////////////////
        sorter.doneAdding();
        progress = new ProgressLogger(logger, 1000000, "written");
        logger.info("Writing out sorted records to final VCF.");

        for (final VariantContext ctx : sorter) {
            out.add(ctx);
            progress.record(ctx.getContig(), ctx.getStart());
        }
        out.close();
        sorter.cleanup();

        return null;
    }
}
