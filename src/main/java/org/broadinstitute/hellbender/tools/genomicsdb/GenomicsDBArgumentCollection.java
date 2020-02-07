package org.broadinstitute.hellbender.tools.genomicsdb;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;

import java.io.Serializable;

public class GenomicsDBArgumentCollection extends GenomicsDBBaseArgumentCollection {
  private static final long serialVersionUID = 1L;

  private static final String MAX_ALTERNATE_ALLELES_LONG_NAME = "genomicsdb-max-alternate-alleles";
  private static final String MAX_GENOTYPE_COUNT_LONG_NAME = "genomicsdb-max-genotype-count";

  @Advanced
  @Argument(fullName=MAX_ALTERNATE_ALLELES_LONG_NAME, doc="Maximum number of alternate alleles to genotype", optional=true)
  public int maxAlternateAlleles = GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;

  @Advanced
  @Argument(fullName=MAX_GENOTYPE_COUNT_LONG_NAME, doc="Maximum number of genotypes to consider at any site", optional=true)
  public int maxGenotypeCount = GenotypeCalculationArgumentCollection.DEFAULT_MAX_GENOTYPE_COUNT;
}