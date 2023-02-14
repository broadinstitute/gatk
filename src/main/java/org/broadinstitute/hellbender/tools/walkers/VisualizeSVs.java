package org.broadinstitute.hellbender.tools.walkers;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.util.SequenceUtil;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import picard.cmdline.programgroups.VariantEvaluationProgramGroup;

import java.io.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Viewer for structural variants in long reads.
 * <p>
 * For each SV described by a VCF, this program:
 * <ul>
 *     <li>alters a few Kb of reference to match the described SV and writes it to a FASTA file</li>
 *     <li>writes the reads overlapping the SV to a FASTQ file</li>
 *     <li>runs minimap2 to align the reads to the altered reference</li>
 *     <li>executes igv to show the alignments</li>
 * </ul>
 * </p>
 * <p>
 *     On the IGV screen, look at the name of the contig being displayed.  It's a description of the
 *     event in the form "contig_position_svtype_svlen".  The first novel adjacency will be placed
 *     at position 1000 on this altered reference fragment.
 * </p>
 * <h3>Requirements</h3>
 * <ul><li>You must have samtools and minimap2 installed and on your PATH.</il></ul>
 * </p>
 * <h3>Input</h3>
 * <ul>
 * <li>A VCF containing structural variants</li>
 * <li>A SAM/BAM/CRAM file with long reads</li>
 * <li>The reference used by the VCF and the reads</li>
 * <li>A region of the VCF for which an IGV screen will be popped for each variant</li>
 * </ul>
 * <p/>
 * <h3>Output</h3>
 * <ul><li>IGV screen popups showing alignments of reads to a structural variation</li></ul>
 * <p/>
 * <h3>Usage example:</h3>
 * <pre>
 * gatk VisualizeSVs \
 *   -V sv.vcf.gz \
 *   -I reads.bam \
 *   -R ref.fa.gz \
 *   -L chr1:16725245
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "For each SV, remap overlapping reads to a revised ref implied by the alt allele, and invoke IGV.",
        oneLineSummary = "Look at how the reads map to a revised reference.",
        programGroup = VariantEvaluationProgramGroup.class
)
public final class VisualizeSVs extends VariantWalker {
    public static final int WINDOW_SIZE = 1000;
    public static final int HUGE_EVENT_SIZE = 10000;

    // ]p]t -- the BND locus, p, is the end of the reference sequence that precedes t
    public static final Pattern BND_PRE_FWD = Pattern.compile("]([!-~]*):([1-9][0-9]*)]([acgtnACGTN]+)");

    // [p[t -- the BND locus, p, is the start of the reference sequence the reverse complement of which precedes t
    public static final Pattern BND_PRE_REV = Pattern.compile("\\[([!-~]*):([1-9][0-9]*)\\[([acgtnACGTN]+)");

    // t[p[ -- the BND locus, p, is the start of the reference sequence that follows t
    public static final Pattern BND_POST_FWD = Pattern.compile("([acgtnACGTN]+)\\[([!-~]*):([1-9][0-9]*)\\[");

    // t]p] -- the BND locus, p, is the end of the reference sequence the reverse complement of which follows t
    public static final Pattern BND_POST_REV = Pattern.compile("([acgtnACGTN]+)]([!-~]*):([1-9][0-9]*)]");

    public static final String scriptText =
            "#!/bin/sh\n" +
            "samtools faidx ref.fa &&\\\n" +
            "    minimap2 -axmap-hifi ref.fa reads.fq | samtools sort -OBAM - > align.bam &&\\\n" +
            "    samtools index align.bam &&\\\n" +
            "    igv -g ref.fa align.bam\n";

    public final File tmpDir = IOUtils.createTempDir("vsv");
    public final File refFile = new File(tmpDir, "ref.fa");
    public final File readsFile = new File(tmpDir, "reads.fq");
    public final File scriptFile = new File(tmpDir, "run.sh");
    public final ProcessBuilder scriptRunner = new ProcessBuilder(scriptFile.getAbsolutePath());

    @Override
    public void onTraversalStart() {
        super.onTraversalStart();
        IOUtils.writeByteArrayToFile(scriptText.getBytes(), scriptFile);
        if ( !scriptFile.setExecutable(true) ) {
            throw new UserException("Can't make minimap2/igv script executable.");
        }
        scriptRunner.directory(tmpDir).inheritIO();
    }

    @Override
    public void apply( final VariantContext variant,
                       final ReadsContext readsContext,
                       final ReferenceContext referenceContext,
                       final FeatureContext featureContext) {
        final String svType = variant.getAttributeAsString(GATKSVVCFConstants.SVTYPE, null);
        final String contig = variant.getContig();
        final int start = variant.getStart();
        final int svLen = variant.getAttributeAsInt(GATKSVVCFConstants.SVLEN, 0);
        final String seqName = contig + "_" + start + "_" + svType + "_" + Math.abs(svLen);

        final int leadStart = Math.max(1, start - WINDOW_SIZE + 1);
        final SimpleInterval leadInterval = new SimpleInterval(contig, leadStart, start);
        final String leadBases = new String(referenceContext.getBases(leadInterval));

        final int end = variant.getEnd();
        final SimpleInterval lagInterval = new SimpleInterval(contig, end + 1, end + WINDOW_SIZE);
        final String lagBases = new String(referenceContext.getBases(lagInterval));

        if ( GATKSVVCFConstants.SYMB_ALT_STRING_DEL.equals(svType) ) {
            final String refSeq = leadBases + lagBases;
            writeRef( seqName, refSeq, false );
        } else if ( GATKSVVCFConstants.SYMB_ALT_STRING_INS.equals(svType) ) {
            final String altBases = variant.getAlternateAllele(0).getBaseString();
            if ( svLen < HUGE_EVENT_SIZE ) {
                final String refSeq = leadBases + altBases.substring(1) + lagBases;
                writeRef(seqName, refSeq, false);
            } else {
                final String bnd1Bases = leadBases + altBases.substring(1, WINDOW_SIZE + 1);
                writeRef(seqName + "_BND1", bnd1Bases, false);
                final int len = altBases.length();
                final String bnd2Bases = altBases.substring(len - WINDOW_SIZE, len) + lagBases;
                writeRef(seqName + "_BND2", bnd2Bases, true);
            }
        } else if ( GATKSVVCFConstants.SYMB_ALT_STRING_DUP.equals(svType) ) {
            final SimpleInterval dupInterval = new SimpleInterval(contig, start + 1, end);
            final String dupBases = new String(referenceContext.getBases(dupInterval));
            if ( svLen < HUGE_EVENT_SIZE ) {
                final String refSeq = leadBases + dupBases + dupBases + lagBases;
                writeRef(seqName, refSeq, false);
            } else {
                final int len = dupBases.length();
                final String refSeq = dupBases.substring(len - WINDOW_SIZE, len) +
                                        dupBases.substring(0, WINDOW_SIZE);
                writeRef(seqName + "_BND", refSeq, false);
            }
        } else if ( GATKSVVCFConstants.SYMB_ALT_STRING_INV.equals(svType) ) {
            final SimpleInterval invInterval = new SimpleInterval(contig, start+1, end);
            final byte[] invBaseCalls = referenceContext.getBases(invInterval);
            SequenceUtil.reverseComplement(invBaseCalls);
            final String invBases = new String(invBaseCalls);
            if ( svLen < HUGE_EVENT_SIZE ) {
                final String refSeq = leadBases + invBases + lagBases;
                writeRef(seqName, refSeq, false);
            } else {
                final String bnd1Bases = leadBases + invBases.substring(0, WINDOW_SIZE);
                writeRef(seqName + "_BND1", bnd1Bases, false);
                final int len = invBases.length();
                final String bnd2Bases = invBases.substring(len - WINDOW_SIZE, len) + lagBases;
                writeRef(seqName + "_BND2", bnd2Bases, true);
            }
        } else if ( GATKSVVCFConstants.BREAKEND_STR.equals(svType) ) {
            final String bndDescription = variant.getAlternateAllele(0).getDisplayString();
            Matcher matcher;
            if ( (matcher = BND_PRE_FWD.matcher(bndDescription)).matches() ) {
                final String tig = matcher.group(1);
                final int bndPosEnd = Integer.parseInt(matcher.group(2));
                final int bndPosBeg = Math.max(1, bndPosEnd - WINDOW_SIZE + 1);
                final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                final String bndBases = new String(referenceContext.getBases(bndInterval));
                final String refSeq = bndBases + matcher.group(3) + lagBases.substring(1);
                writeRef( seqName, refSeq, false );
            } else if ( (matcher = BND_PRE_REV.matcher(bndDescription)).matches() ) {
                final String tig = matcher.group(1);
                final int bndPosBeg = Integer.parseInt(matcher.group(2));
                final int bndPosEnd = bndPosBeg + WINDOW_SIZE - 1;
                final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                final byte[] bndBaseCalls = referenceContext.getBases(bndInterval);
                SequenceUtil.reverseComplement(bndBaseCalls);
                final String refReplace = matcher.group(3);
                final String refSeq = new String(bndBaseCalls) + refReplace + lagBases.substring(1);
                writeRef( seqName, refSeq, false );
            } else if ( (matcher = BND_POST_FWD.matcher(bndDescription)).matches() ) {
                final String refReplace = matcher.group(1);
                final String tig = matcher.group(2);
                final int bndPosBeg = Integer.parseInt(matcher.group(3));
                final int bndPosEnd = bndPosBeg + WINDOW_SIZE - 1;
                final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                final String bndBases = new String(referenceContext.getBases(bndInterval));
                final String refSeq =
                        leadBases.substring(0, leadBases.length() - 1) + refReplace + bndBases;
                writeRef( seqName, refSeq, false );
            } else if ( (matcher = BND_POST_REV.matcher(bndDescription)).matches() ) {
                final String refReplace = matcher.group(1);
                final String tig = matcher.group(2);
                final int bndPosEnd = Integer.parseInt(matcher.group(3));
                final int bndPosBeg = Math.max(1, bndPosEnd - WINDOW_SIZE + 1);
                final SimpleInterval bndInterval = new SimpleInterval(tig, bndPosBeg, bndPosEnd);
                final byte[] bndBaseCalls = referenceContext.getBases(bndInterval);
                SequenceUtil.reverseComplement(bndBaseCalls);
                final String refSeq = leadBases.substring(0, leadBases.length() - 1) +
                        refReplace + new String(bndBaseCalls);
                writeRef( seqName, refSeq, false );
            } else {
                logger.warn("Can't interpret the BND description: " + bndDescription);
                return;
            }
        } else if ( svType != null ) {
            logger.warn("Don't know how to handle SVTYPE=" + svType);
            return;
        }
        writeReads( readsContext );
        alignAndDisplay();
    }

    private void alignAndDisplay() {
        try {
            final Process process = scriptRunner.start();
            final int exitValue = process.waitFor();
            if ( exitValue != 0 ) {
                throw new UserException("Script that runs minimap2 and igv returned " + exitValue);
            }
        } catch ( final IOException | InterruptedException ex ) {
            throw new UserException("Script that runs minimap2 and igv failed.", ex);
        }
    }

    private void writeRef( final String seqName, final String bases, final boolean append ) {
        try ( final BufferedWriter writer =
                  new BufferedWriter(new OutputStreamWriter(new FileOutputStream(refFile, append))) ) {
            writer.write(">");
            writer.write(seqName);
            writer.newLine();
            writeLines(writer, bases);
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write altered reference.", ioe);
        }
    }

    private void writeReads( final ReadsContext reads ) {
        try ( final BufferedWriter writer =
                  new BufferedWriter(new OutputStreamWriter(new FileOutputStream(readsFile))) ) {
            for ( final GATKRead read : reads ) {
                writer.write("@");
                writer.write(read.getName());
                writer.newLine();
                writeLines(writer, read.getBasesString());
                writer.write("+");
                writer.newLine();
                writeLines(writer, SAMUtils.phredToFastq(read.getBaseQualitiesNoCopy()));
            }
        } catch ( final IOException ioe ) {
            throw new UserException("Can't write reads.", ioe);
        }
    }

    private void writeLines( final BufferedWriter writer, final String seq ) throws IOException {
        final int nBases = seq.length();
        int beg = 0;
        while ( beg < nBases ) {
            final int end = Math.min(nBases, beg + 80);
            writer.write(seq.substring(beg, end));
            writer.newLine();
            beg = end;
        }
    }
}
