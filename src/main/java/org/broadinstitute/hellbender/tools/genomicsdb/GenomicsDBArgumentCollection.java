package org.broadinstitute.hellbender.tools.genomicsdb;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class GenomicsDBArgumentCollection implements Serializable {
  private static final long serialVersionUID = 1L;
  public static final String USE_BCF_CODEC_LONG_NAME = "genomicsdb-use-bcf-codec";
  public static final String SHARED_POSIXFS_OPTIMIZATIONS = "genomicsdb-shared-posixfs-optimizations";

  public static final String CALL_GENOTYPES_LONG_NAME = "call-genotypes";
  private static final boolean DEFAULT_CALL_GENOTYPES = false;
  private static final boolean DEFAULT_USE_BCF_CODEC = false;
  private static final boolean DEFAULT_SHARED_POSIXFS_OPTIMIZATIONS = false;

  /**
   * Not full-fledged arguments for now.
   */
  public int maxDiploidAltAllelesThatCanBeGenotyped = GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;

  @Argument(fullName = CALL_GENOTYPES_LONG_NAME, doc = "Output called genotypes in final VCF (otherwise no-call)", optional = true)
  public boolean callGenotypes = DEFAULT_CALL_GENOTYPES;

  /**
   * Currently there is no support for 64-bit fields in BCF2Codec. The VCFCodec allows for 64-bit
   * width positions and INFO fields and for computed annotation sizes to exceed the 32-bit
   * integer space while encoding/decoding with GenomicsDB. Use the BCF2Codec option if and
   * only if performance is an issue.
   */
  @Advanced
  @Argument(
      fullName = USE_BCF_CODEC_LONG_NAME,
      doc =
          "Use BCF Codec Streaming for data from GenomicsDB instead of the default VCFCodec. BCFCodec performs slightly better but currently does not support "
              + "64-bit width positions and INFO fields and for computed annotation sizes to exceed 32-bit integer space.",
      optional = true)
  public boolean useBCFCodec = DEFAULT_USE_BCF_CODEC;

  @Argument(fullName = SHARED_POSIXFS_OPTIMIZATIONS,
          doc = "Allow for optimizations to improve the usability and performance for shared Posix Filesystems(e.g. NFS, Lustre). " +
                  "If set, file level locking is disabled and file system writes are minimized.",
          optional = true)
  public boolean sharedPosixFSOptimizations = DEFAULT_SHARED_POSIXFS_OPTIMIZATIONS;
}
