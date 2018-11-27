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
import org.broadinstitute.gatk.utils.commandline.ArgumentCollection;
import org.broadinstitute.gatk.utils.commandline.Input;
import org.broadinstitute.gatk.utils.commandline.RodBinding;
import org.broadinstitute.gatk.engine.CommandLineGATK;
import org.broadinstitute.gatk.engine.arguments.StandardVariantContextInputArgumentCollection;
import org.broadinstitute.gatk.utils.contexts.AlignmentContext;
import org.broadinstitute.gatk.utils.contexts.ReferenceContext;
import org.broadinstitute.gatk.utils.refdata.RefMetaDataTracker;
import org.broadinstitute.gatk.engine.walkers.*;
import org.broadinstitute.gatk.utils.BaseUtils;
import org.broadinstitute.gatk.utils.GenomeLoc;
import org.broadinstitute.gatk.engine.SampleUtils;
import org.broadinstitute.gatk.utils.collections.Pair;
import org.broadinstitute.gatk.utils.exceptions.UserException;
import org.broadinstitute.gatk.utils.help.DocumentedGATKFeature;
import org.broadinstitute.gatk.utils.help.HelpConstants;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.gatk.utils.variant.GATKVCFConstants;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.Optional;


/**
 * Generate an alternative reference sequence over the specified interval
 *
 * <p>Given a variant callset, this tool replaces the reference bases at variation sites with the bases supplied in the
 * corresponding callset records. Additionally, it allows for one or more "snpmask" VCFs to set overlapping bases to 'N'.</p>
 *
 * <p>The output format can be partially controlled using the provided command-line arguments.
 * Specify intervals with the usual -L argument to output only the reference bases within your intervals.
 * Overlapping intervals are automatically merged; reference bases for each disjoint interval will be output as a
 * separate fasta sequence (named numerically in order).</p>
 *
 * <h3>Caveats</h3>
 * <ul>
 *     <li>If there are multiple variants that start at a site, it chooses one of them randomly.</li>
 *     <li>When there are overlapping indels (but with different start positions) only the first will be chosen.</li>
 *     <li>This tool works only for SNPs and for simple indels (but not for things like complex substitutions).</li>
 * </ul>

 * <h3>Input</h3>
 * <p>
 * The reference, requested intervals, and any number of variant ROD files.
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A FASTA file representing the requested intervals.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 * java -jar GenomeAnalysisTK.jar \
 *   -T FastaAlternateReferenceMaker \
 *   -R reference.fasta \
 *   -o output.fasta \
 *   -L input.intervals \
 *   -V input.vcf \
 *   [--snpmask mask.vcf]
 * </pre>
 *
 */
@DocumentedGATKFeature( groupName = HelpConstants.DOCS_CAT_REFUTILS, extraDocs = {CommandLineGATK.class} )
@Reference(window=@Window(start=-1,stop=50))
@Requires(value={DataSource.REFERENCE})
public class FastaAlternateReferenceMaker extends FastaReferenceMaker {

  /**
   * Variants from this input file are used by this tool to construct an alternate reference.
   */
  @ArgumentCollection
  protected StandardVariantContextInputArgumentCollection variantCollection = new StandardVariantContextInputArgumentCollection();

  /**
   * SNPs from this file are used as a mask (inserting N's in the sequence) when constructing the alternate reference
   */
  @Input(fullName="snpmask", shortName = "snpmask", doc="SNP mask VCF file", required=false)
  protected RodBinding<VariantContext> snpmask;

  /**
   * Gives priority to a SNP mask over an input VCF for a site. Only has an effect if the --snpmask argument is used.
   */
  @Argument(fullName="snpmaskPriority", shortName = "snpmaskPriority", doc="SNP mask priority", required=false)
  protected Boolean snpmaskPriority = false;

  /**
   * This option will generate an error if the specified sample does not exist in the VCF.
   * Non-diploid (or non-called) genotypes are ignored.
   */
  @Argument(fullName="use_IUPAC_sample", shortName="IUPAC", doc = "If specified, heterozygous SNP sites will be output using IUPAC ambiguity codes given the genotypes for this sample", required=false)
  private String iupacSample = null;

  private int deletionBasesRemaining = 0;

  private static final String EMPTY_BASE = " ";

  @Override
  public void initialize() {
    super.initialize();
    if ( iupacSample != null ) {
      final List<String> rodName = Arrays.asList(variantCollection.variants.getName());
      final Set<String> samples = SampleUtils.getUniqueSamplesFromRods(getToolkit(), rodName);
      if ( !samples.contains(iupacSample) )
        throw new UserException.BadInput("the IUPAC sample specified is not present in the provided VCF file");
    }
  }

  @Override
  public Pair<GenomeLoc, String> map(final RefMetaDataTracker tracker, final ReferenceContext ref, final AlignmentContext context) {

    if (deletionBasesRemaining > 0) {
      deletionBasesRemaining--;
      return new Pair<>(context.getLocation(), "");
    }

    final String refBase = String.valueOf((char)ref.getBase());

    // If we have a mask at this site, use it
    if ( snpmaskPriority ){
      final Pair<GenomeLoc, String> mask = maskSnp(tracker, context);
      if ( mask != null )
        return mask;
    }

    // Check to see if we have a called snp
    for ( final VariantContext vc : tracker.getValues(variantCollection.variants, ref.getLocus()) ) {
      if ( vc.isFiltered() )
        continue;

      if ( vc.isSimpleDeletion()) {
        deletionBasesRemaining = vc.getReference().length() - 1;
        // delete the next n bases, not this one
        return new Pair<>(context.getLocation(), refBase);
      } else if ( vc.isSimpleInsertion() || vc.isSNP() ) {
        // Get the first alt allele that is not a spanning deletion. If none present, use the empty allele
        final Optional<Allele> optionalAllele = getFirstNonSpanDelAltAllele(vc.getAlternateAlleles());
        final Allele allele = optionalAllele.isPresent() ? optionalAllele.get() : Allele.create(EMPTY_BASE, false);
        if ( vc.isSimpleInsertion() ) {
          return new Pair<>(context.getLocation(), allele.toString());
        } else {
          final String base = (iupacSample != null) ? getIUPACbase(vc.getGenotype(iupacSample), refBase) : allele.toString();
          return new Pair<>(context.getLocation(), base);
        }
      }
    }

    if ( !snpmaskPriority ){
      final Pair<GenomeLoc, String> mask = maskSnp(tracker, context);
      if ( mask != null )
        return mask;
    }

    // if we got here then we're just ref
    return new Pair<>(context.getLocation(), refBase);
  }

  /**
   * Get the first non spanning deletion (* or <*:DEL>) alt allele
   * @param altAlleles the alternate alleles
   * @return the first non spanning deletion allele or null
   */
  private Optional<Allele> getFirstNonSpanDelAltAllele( final List<Allele> altAlleles ) {
    for (final Allele allele : altAlleles) {
      if (!allele.equals(Allele.SPAN_DEL) && !allele.equals(GATKVCFConstants.SPANNING_DELETION_SYMBOLIC_ALLELE_DEPRECATED)) {
        return Optional.of(allele);
      }
    }

    return Optional.empty();
  }

  /**
   * Mask a SNP (inserting N's in the sequence)
   *
   * @param tracker the Reference Metadata available at a particular site in the genome
   * @param context the locus context data
   * @return mask at the locus or null if no SNP at that locus
   */
  private Pair<GenomeLoc, String> maskSnp(final RefMetaDataTracker tracker, final AlignmentContext context){
    for (final VariantContext vc : tracker.getValues(snpmask)) {
      if (vc.isSNP()) {
        return new Pair<>(context.getLocation(), "N");
      }
    }

    return null;
  }

  /**
   * Returns the IUPAC encoding for the given genotype or the reference base if not possible
   *
   * @param genotype  the genotype to encode
   * @param ref       the reference base
   * @return non-null, non-empty String of bases
   */
  private String getIUPACbase(final Genotype genotype, final String ref) {
    if ( genotype == null )
      throw new IllegalStateException("The genotype is null for sample " + iupacSample);

    // If have a spanning deletion, if both alleles are spanning deletions, use the empty allele. Otherwise, use the allele is not a
    // spanning deletion.
    if ( genotype.getAlleles().contains(Allele.SPAN_DEL) ) {
      if ( genotype.isHomVar() ) {
        return EMPTY_BASE;
      } else {
        return genotype.getAllele(0).equals(Allele.SPAN_DEL) ? genotype.getAllele(1).getBaseString() : genotype.getAllele(0).getBaseString();
      }
    }

    if ( !genotype.isHet() )
      return genotype.isHom() ? genotype.getAllele(0).getBaseString() : ref;

    final byte allele1 = genotype.getAllele(0).getBases()[0];
    final byte allele2 = genotype.getAllele(1).getBases()[0];
    return new String(new byte[] {BaseUtils.basesToIUPAC(allele1, allele2)});
  }
}