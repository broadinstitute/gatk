package org.broadinstitute.hellbender.tools.genomicsdb;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection;

public class GenomicsDBArgumentCollection extends GenomicsDBBaseArgumentCollection {
  private static final long serialVersionUID = 1L;

  private static final String MAX_GENOTYPE_COUNT_LONG_NAME = "genomicsdb-max-genotype-count";

  @Advanced
  @Argument(fullName=MAX_GENOTYPE_COUNT_LONG_NAME, doc="Maximum number of genotypes to consider at any site", optional=true)
  public int maxGenotypeCount = GenotypeCalculationArgumentCollection.DEFAULT_MAX_GENOTYPE_COUNT;
}