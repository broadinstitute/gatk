package org.broadinstitute.hellbender.tools.genomicsdb;

import htsjdk.variant.variantcontext.GenotypeLikelihoods;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;

public class GenomicsDBArgumentCollection implements Serializable {
  private static final long serialVersionUID = 1L;
  public static final String USE_BCF_CODEC_LONG_NAME = "genomicsdb-use-bcf-codec";
  public static final String SHARED_POSIXFS_OPTIMIZATIONS = "genomicsdb-shared-posixfs-optimizations";
  public static final String USE_GCS_HDFS_CONNECTOR = "genomicsdb-use-gcs-hdfs-connector";

  public static final String CALL_GENOTYPES_LONG_NAME = "call-genotypes";
  public static final String MAX_ALTS_LONG_NAME = "genomicsdb-max-alternate-alleles";
  private static final boolean DEFAULT_CALL_GENOTYPES = false;
  private static final boolean DEFAULT_USE_BCF_CODEC = false;
  private static final boolean DEFAULT_SHARED_POSIXFS_OPTIMIZATIONS = false;
  private static final boolean DEFAULT_USE_GCS_HDFS_CONNECTOR = false;

  /**
   * Maximum number of alternate alleles that will report likelihoods after being combined on reading from GenomicsDB (including <NON_REF>)
   * Must be at least one greater than the maximum number of alternate alleles for genotyping.
   * A typical value is 3 more than the --max-alternate-alleles value that's used by GenotypeGVCFs and larger differences
   * result in more robustness to PCR-related indel errors.
   * Note that GenotypeGVCFs will drop highly multi-allelic sites that are missing likelihoods.
   *
   * See also {@link org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeCalculationArgumentCollection#MAX_ALTERNATE_ALLELES_LONG_NAME}
   */
  @Argument(fullName = MAX_ALTS_LONG_NAME, doc = "Maximum number of alternate alleles that will be combined on reading from GenomicsDB")
  public int maxDiploidAltAllelesThatCanBeGenotyped = GenotypeLikelihoods.MAX_DIPLOID_ALT_ALLELES_THAT_CAN_BE_GENOTYPED;

    /**
     * Output called genotypes in the final VCF (otherwise no-call)
     */
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

  @Argument(fullName = USE_GCS_HDFS_CONNECTOR,
          doc = "Use the GCS HDFS Connector instead of the native GCS SDK client with gs:// URLs.",
          optional = true)
  public boolean useGcsHdfsConnector = DEFAULT_USE_GCS_HDFS_CONNECTOR;
}
