package org.broadinstitute.hellbender.tools.genomicsdb;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class GenomicsDBArgumentCollection implements Serializable {
  private static final long serialVersionUID = 1L;
  public static final String USE_VCF_CODEC_LONG_NAME = "genomicsdb-use-vcf-codec";

  private static final boolean DEFAULT_CALL_GENOTYPES = false;
  private static final boolean DEFAULT_USE_VCF_CODEC = false;

  /**
   * Not full-fledged arguments for now.
   */
  public boolean callGenotypes = DEFAULT_CALL_GENOTYPES;
  public int maxDiploidAltAllelesThatCanBeGenotyped = GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;


  /**
   * Currently there is no support for 64-bit fields in BCFCodec. Specifying this option will allow
   * for 64-bit width positions and INFO fields and for computed annotation sizes to exceed the 32-bit
   * integer space while encoding/decoding with GenomicsDB.
   */
  @Advanced
  @Argument(
      fullName = USE_VCF_CODEC_LONG_NAME,
      doc = "Use VCF Codec Streaming for data from GenomicsDB instead of the default BCF",
      optional = true)
  public boolean useVCFCodec = DEFAULT_USE_VCF_CODEC;
}