/*
* Copyright 2012-2015 Broad Institute, Inc.
* 
* Permission is hereby granted, free of charge, to any person
* obtaining a copy of this software and associated documentation
* files (the "Software"), to deal in the Software without
* restriction, including without limitation the rights to use,
* copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the
* Software is furnished to do so, subject to the following
* conditions:
* 
* The above copyright notice and this permission notice shall be
* included in all copies or substantial portions of the Software.
* 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
* OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
* NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
* HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
* WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
* THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

package org.broadinstitute.gatk.tools.walkers.qc;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.DataSource;
import org.broadinstitute.gatk.engine.walkers.LocusWalker;
import org.broadinstitute.gatk.engine.walkers.Requires;
import org.broadinstitute.gatk.engine.walkers.TreeReducible;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.codecs.sampileup.SAMPileupFeature;
import org.broadinstitute.gatk.utils.exceptions.ReviewedGATKException;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import org.broadinstitute.gatk.utils.pileup.ReadBackedPileup;

import java.io.PrintStream;
import java.util.Arrays;

/**
 * Compare GATK's internal pileup to a reference Samtools pileup
 *
 * <p>At every locus in the input set, compares the pileup data (reference base, aligned base from
 * each overlapping read, and quality score) generated internally by GATK to a reference pileup data generated
 * by Samtools. Note that the pileup program has been replaced in Samtools by mpileup, which produces a slightly
 * different output format by default.
 * </p>
 *
 * <h3>Format</h3>
 * <p>There are two versions of the original pileup format: the current 6-column format produced by Samtools, and the old
 * 10-column "consensus" format which could be obtained by using the -c argument, now deprecated.</p>
 * <h4>Simple pileup: 6-column format</h4>
 * <p>
 *     Each line consists of chromosome, 1-based coordinate, reference base, the
 *     number of reads covering the site, read bases and base qualities. At the
 *     read base column, a dot stands for a match to the reference base on the
 *     forward strand, a comma for a match on the reverse strand, `ACGTN' for a mismatch
 *     on the forward strand and `acgtn' for a mismatch on the reverse strand.
 *     A pattern `\+[0-9]+[ACGTNacgtn]+' indicates there is an insertion between
 *     this reference position and the next reference position. The length of the
 *     insertion is given by the integer in the pattern, followed by the inserted sequence.
 * </p>
 * <pre>
 *     seq1 272 T 24  ,.$.....,,.,.,...,,,.,..^+. <<<+;<<<<<<<<<<<=<;<;7<&
 *     seq1 273 T 23  ,.....,,.,.,...,,,.,..A <<<;<<<<<<<<<3<=<<<;<<+
 *     seq1 274 T 23  ,.$....,,.,.,...,,,.,...    7<7;<;<<<<<<<<<=<;<;<<6
 *     seq1 275 A 23  ,$....,,.,.,...,,,.,...^l.  <+;9*<<<<<<<<<=<<:;<<<<
 *     seq1 276 G 22  ...T,,.,.,...,,,.,....  33;+<<7=7<<7<&<<1;<<6<
 *     seq1 277 T 22  ....,,.,.,.C.,,,.,..G.  +7<;<<<<<<<&<=<<:;<<&<
 *     seq1 278 G 23  ....,,.,.,...,,,.,....^k.   %38*<<;<7<<7<=<<<;<<<<<
 *     seq1 279 C 23  A..T,,.,.,...,,,.,..... ;75&<<<<<<<<<=<<<9<<:<<
 * </pre>
 * <p>
 *     See the <a href="http://samtools.sourceforge.net/pileup.shtml">Pileup format documentation</a> for more details.
 * </p>
 *
 * <h4>Consensus pileup: 10/13-column format</h4>
 * <p>The "consensus" or extended pileup consists of the following:
 *  <ul>
 *      <li>original 6 columns as described above</li>
 *      <li>4 extra columns representing consensus values (consensus base, consensus quality, variant quality and maximum mapping quality of the
 * reads covering the sites) for all sites, inserted before the bases and quality strings</li>
 *      <li>3 extra columns indicating counts of reads supporting indels (just for indel sites)</li>
 *  </ul>
 * </p>
 * <h4>Example of consensus pileup for SNP or non-variant sites</h4>
 * <pre>
 *     seq1  60  T  T  66  0  99  13  ...........^~.^~.   9<<55<;<<<<<<
 *     seq1  61  G  G  72  0  99  15  .............^~.^y. (;975&;<<<<<<<<
 *     seq1  62  T  T  72  0  99  15  .$..............    <;;,55;<<<<<<<<
 *     seq1  63  G  G  72  0  99  15  .$.............^~.  4;2;<7:+<<<<<<<
 *     seq1  64  G  G  69  0  99  14  ..............  9+5<;;;<<<<<<<
 *     seq1  65  A  A  69  0  99  14  .$............. <5-2<;;<<<<<<;
 *     seq1  66  C  C  66  0  99  13  .............   &*<;;<<<<<<8<
 *     seq1  67  C  C  69  0  99  14  .............^~.    ,75<.4<<<<<-<<
 *     seq1  68  C  C  69  0  99  14  ..............  576<;7<<<<<8<< *
 * </pre>
 *
 * <h4>Example of consensus pileup for indels</h4>
 * <pre>
 *     Escherichia_coli_K12	3995037	*	*\/*	430	0	37	144	*	+A	143	1	0
 *     Escherichia_coli_K12	3995279	*	*\/*	202	0	36	68	*	+A	67	1	0
 *     Escherichia_coli_K12	3995281	*	*\/*	239	0	36	67	*	-CG	66	1	0
 * </pre>
 * <p>
 *     See <a href="http://samtools.sourceforge.net/cns0.shtml/">Consensus pileup format (deprecated)</a> for more details.
 * </p>
 *
 * <h3>Input</h3>
 * <p>A BAM file containing your aligned sequence data and a pileup file generated by Samtools covering the region you
 * want to examine.</p>
 *
 * <h3>Output</h3>
 * <p>A text file listing mismatches between the input pileup and the GATK's internal pileup. If there are no mismatches, the output file is empty.</p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T CheckPileup \
 *   -R reference.fasta \
 *   -I your_data.bam \
 *   --pileup:SAMPileup pileup_file.txt \
 *   -L chr1:257-275 \
 *   -o output_file_name
 * </pre>
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_QC, extraDocs = {CommandLineGATK.class} )
@Requires(value={DataSource.READS,DataSource.REFERENCE})
public class CheckPileup extends LocusWalker<Integer, CheckPileupStats> implements TreeReducible<CheckPileupStats> {
    /**
     * This is the existing pileup against which we'll compare GATK's internal pileup at each genome position in the desired interval.
     */
    @Input(fullName = "pileup", shortName = "pileup", doc="Pileup generated by Samtools", required = true)
    RodBinding<SAMPileupFeature> pileup;

    @Output
    private PrintStream out;
    /**
     * By default the program will quit if it encounters an error (such as missing truth data for a given position).
     * Use this flag to override the default behavior; the program will then simply print an error message and move on
     * to the next position.
     */
    @Argument(fullName="continue_after_error",doc="Continue after encountering an error",required=false)
    public boolean CONTINUE_AFTER_AN_ERROR = false;

    public Integer map(RefMetaDataTracker tracker, ReferenceContext ref, AlignmentContext context) {
        ReadBackedPileup pileup = context.getBasePileup();
        SAMPileupFeature truePileup = getTruePileup( tracker );

        if ( truePileup == null ) {
            out.printf("No truth pileup data available at %s%n", pileup.getPileupString(ref.getBaseAsChar()));
            if ( ! CONTINUE_AFTER_AN_ERROR ) {
                throw new UserException.BadInput(String.format("No pileup data available at %s given GATK's output of %s -- this walker requires samtools pileup data over all bases",
                        context.getLocation(), new String(pileup.getBases())));
            }
        } else {
            String pileupDiff = pileupDiff(pileup, truePileup, true);
            if ( pileupDiff != null ) {
                out.printf("%s vs. %s%n", pileup.getPileupString(ref.getBaseAsChar()), truePileup.getPileupString());
                if ( ! CONTINUE_AFTER_AN_ERROR ) {
                    throw new UserException.BadInput(String.format("The input pileup doesn't match the GATK's internal pileup: %s", pileupDiff));
                }
            }
        }

        return pileup.getNumberOfElements();
    }

    private static String maybeSorted( final String x, boolean sortMe )
    {
        if ( sortMe ) {
            byte[] bytes = x.getBytes();
            Arrays.sort(bytes);
            return new String(bytes);
        }
        else
            return x;
    }

    public String pileupDiff(final ReadBackedPileup a, final SAMPileupFeature b, boolean orderDependent)
    {
        if ( a.getNumberOfElements() != b.size() )
            return "Sizes not equal";
        GenomeLoc featureLocation = getToolkit().getGenomeLocParser().createGenomeLoc(b.getChr(),b.getStart(),b.getEnd());
        if ( a.getLocation().compareTo(featureLocation) != 0 )
            return "Locations not equal";

        String aBases = maybeSorted(new String(a.getBases()), ! orderDependent );
        String bBases = maybeSorted(b.getBasesAsString(), ! orderDependent );
        if ( ! aBases.toUpperCase().equals(bBases.toUpperCase()) )
            return "Bases not equal";

        String aQuals = maybeSorted(new String(a.getQuals()), ! orderDependent );
        String bQuals = maybeSorted(new String(b.getQuals()), ! orderDependent );
        if ( ! aQuals.equals(bQuals) )
            return "Quals not equal";

        return null;
    }

    // Given result of map function
    public CheckPileupStats reduceInit() { return new CheckPileupStats(); }
    public CheckPileupStats reduce(Integer value, CheckPileupStats sum) {
        sum.nLoci++;
        sum.nBases += value;
        return sum;
    }

    public CheckPileupStats treeReduce( CheckPileupStats lhs, CheckPileupStats rhs ) {
        CheckPileupStats combined = new CheckPileupStats();
        combined.nLoci = lhs.nLoci + rhs.nLoci;
        combined.nBases = lhs.nBases + rhs.nBases;
        return combined;
    }

    /**
     * Extracts the true pileup data from the given rodSAMPileup.  Note that this implementation
     * assumes that the genotype will only be point or indel.
     * @param tracker ROD tracker from which to extract pileup data.
     * @return True pileup data.
     */
    private SAMPileupFeature getTruePileup( RefMetaDataTracker tracker ) {
        SAMPileupFeature pileupArg = tracker.getFirstValue(pileup);

        if( pileupArg == null)
            return null;

        if( pileupArg.hasPointGenotype() )
            return pileupArg.getPointGenotype();
        else if( pileupArg.hasIndelGenotype() )
            return pileupArg.getIndelGenotype();
        else
            throw new ReviewedGATKException("Unsupported pileup type: " + pileupArg);
    }
}

class CheckPileupStats {
    public long nLoci = 0;
    public long nBases = 0;

    public CheckPileupStats() {
    }

    public String toString() {
        return String.format("Validated %d sites covered by %d bases%n", nLoci, nBases);
    }
}