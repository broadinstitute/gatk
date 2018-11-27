/*
 * Copyright 2012-2016 Broad Institute, Inc.
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

package org.broadinstitute.hellbender.tools.walkers.fasta;

import org.broadinstitute.gatk.utils.commandline.Argument;
import org.broadinstitute.gatk.utils.commandline.Output;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.RefWalker;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;

import java.io.PrintStream;

/**
 * Create a subset of a FASTA reference sequence
 *
 * <p>This tool creates a new reference in FASTA format consisting of only those positions or intervals
 * provided in the input data set. The output format can be partially controlled using the provided command-line
 * arguments. Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).</p>
 *
 * <h3>Input</h3>
 * <p>
 * The reference and requested intervals.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A fasta file representing the requested intervals. Each interval has a description line starting with a greater-than (">") symbol followed by sequence data.
 * The description begins with the contig name followed by the beginning position on the contig.
 * <pre>
 * For example, the fasta file for contig 1 and intervals 1:3-1:4 and 1:6-1:9
 * >1 1:3
 * AT
 * >1 1:6
 * GGGG
 * </pre>
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T FastaReferenceMaker \
 *   -R reference.fasta \
 *   -o output.fasta \
 *   -L input.intervals
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_REFUTILS, extraDocs = {CommandLineGATK.class} )
public class FastaReferenceMaker extends RefWalker<Pair<GenomeLoc, String>, GenomeLoc> {

  @Output PrintStream out;

  @Argument(fullName="lineWidth", shortName="lw", doc="Maximum length of sequence to write per line", required=false)
  public int fastaLineWidth=60;

  /**
   *  Please note that when using this argument adjacent intervals will automatically be merged.
   */
  @Argument(fullName="rawOnelineSeq", shortName="raw", doc="Print sequences with no FASTA header lines, one line per interval (i.e. lineWidth = infinity)", required=false)
  public boolean fastaRawSeqs=false;

  protected FastaSequence fasta;

  public void initialize() {
    if (fastaRawSeqs) fastaLineWidth = Integer.MAX_VALUE;
    fasta = new FastaSequence(out, fastaLineWidth, fastaRawSeqs);
  }

  public Pair<GenomeLoc, String> map(RefMetaDataTracker rodData, ReferenceContext ref, AlignmentContext context) {
    return new Pair<GenomeLoc, String>(context.getLocation(), String.valueOf((char)ref.getBase()));
  }

  public GenomeLoc reduceInit() {
    return null;
  }

  public GenomeLoc reduce(Pair<GenomeLoc, String> value, GenomeLoc sum) {
    if ( value == null )
      return sum;

    // if there is no interval to the left, then this is the first one
    if ( sum == null ) {
      sum = value.first;
      fasta.setName(fasta.getName() + " " + sum.toString());
      fasta.append(value.second);
    }
    // if the intervals are not contiguous, print out the leftmost one and start a new one
    // (end of contig or new interval)
    else if ( value.first.getStart() != sum.getStop() + 1 || ! value.first.getContig().equals(sum.getContig()) ) {
      fasta.flush();
      sum = value.first;
      fasta.setName(fasta.getName() + " " + sum.toString());
      fasta.append(value.second);
    }
    // otherwise, merge them
    else {
      sum = sum.setStop(sum, value.first.getStop());
      fasta.append(value.second);
    }
    return sum;
  }

  public void onTraversalDone(GenomeLoc sum) {
    fasta.flush();
  }
}