package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils.RunningAverage;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodec;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodecFinder;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * <p>Merges locus-sorted LocusDepth evidence files, and calculates the bi-allelic frequency (baf)
 * for each sample and site, and writes these values as a BafEvidence output file.</p>
 * <p>Samples at sites that have too few depth counts (controlled by --min-total-depth), and samples
 * failing a Pearson's chi square test for goodness of fit to a bi-allelic model are excluded. (I.e.,
 * only samples that are judged to be heterozygous at a site are included.)</p>
 * <p>Across all evaluable samples, sites where the samples have bafs with too much variance
 * (controlled by --max-std) are excluded.  The median baf value at each site is adjusted to be
 * exactly 0.5.</p>
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         One or more LocusDepth evidence files, or a file containing a list of evidence files, one per line.
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         An output file containing merged BafEvidence from the inputs.
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk LocusDepthtoBAF \
 *       -F file1.ld.txt.gz [-F file2.ld.txt.gz ...] \
 *       -O merged.baf.bci
 * </pre>
 *
 * @author Ted Sharpe &lt;tsharpe@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Convert LocusDepth to BafEvidence",
        oneLineSummary = "Convert LocusDepth to BafEvidence",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@ExperimentalFeature
public class LocusDepthtoBAF extends MultiFeatureWalker<LocusDepth> {
    public static final String LOCUS_DEPTH_FILE_NAME = "locus-depth";
    public static final String BAF_SITES_VCF_LONG_NAME = "baf-sites-vcf";
    public static final String SAMPLE_NAMES_NAME = "sample-names";
    public static final String BAF_EVIDENCE_FILE_NAME = "baf-evidence-output";
    public static final String MIN_TOTAL_DEPTH_NAME = "min-total-depth";
    public static final String MAX_STD_NAME = "max-std";
    public static final String MIN_HET_PROBABILITY = "min-het-probability";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";
    public static ChiSquaredDistribution chiSqDist = new ChiSquaredDistribution(null, 1.);

    @Argument(
            doc = "Locus depth files to process",
            fullName = LOCUS_DEPTH_FILE_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME )
    private List<FeatureInput<LocusDepth>> locusDepthFiles;

    @Argument(fullName = BAF_SITES_VCF_LONG_NAME,
            doc = "Input VCF of SNPs marking loci for allele counts")
    public GATKPath bafSitesFile;

    @Argument(fullName = SAMPLE_NAMES_NAME,
            doc = "Sample names. Necessary only when writing a *.baf.bci output file.", optional = true)
    private List<String> sampleNames;

    @Argument(
            doc = "BAF output file",
            fullName = BAF_EVIDENCE_FILE_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME )
    private GATKPath bafOutputFile;

    @Argument(
            doc = "Maximum standard deviation across bafs at a locus (otherwise locus will be ignored)",
            fullName = MAX_STD_NAME, optional = true )
    private double maxStdDev = .2;

    @Argument(
            doc = "Minimum total depth at a locus (otherwise locus will be ignored)",
            fullName = MIN_TOTAL_DEPTH_NAME, optional = true )
    private int minTotalDepth = 10;

    @Argument(
            doc = "Minimum probability of site being biallelic and heterozygous",
            fullName = MIN_HET_PROBABILITY, optional = true )
    private double minHetProbability = .5;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true )
    private int compressionLevel = 4;

    private FeatureDataSource<VariantContext> bafSitesSource;
    private Iterator<VariantContext> bafSitesIterator;
    private FeatureSink<BafEvidence> writer;
    private final List<LocusDepth> sameLocusBuffer = new ArrayList<>();

    @Override
    @SuppressWarnings("unchecked")
    public void onTraversalStart() {
        super.onTraversalStart();
        bafSitesSource = new FeatureDataSource<>(bafSitesFile.toPath().toString());
        bafSitesIterator = bafSitesSource.iterator();
        final SAMSequenceDictionary dict = getDictionary();
        dict.assertSameDictionary(bafSitesSource.getSequenceDictionary());
        final FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>> codec =
                FeatureOutputCodecFinder.find(bafOutputFile);
        if ( !codec.getFeatureType().isAssignableFrom(BafEvidence.class) ) {
            throw new UserException("We're intending to write BafEvidence, but the feature type " +
                    "associated with the output file expects features of type " +
                    codec.getFeatureType().getSimpleName());
        }
        writer = (FeatureSink<BafEvidence>)codec.makeSink(bafOutputFile, dict, sampleNames,
                                                            compressionLevel);
    }

    @Override
    public void apply( final LocusDepth locusDepth, final Object header,
                       final ReadsContext reads, final ReferenceContext ref ) {
        if ( !sameLocus(locusDepth) ) {
            processBuffer();
        }
        sameLocusBuffer.add(locusDepth);
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        if ( !sameLocusBuffer.isEmpty() ) {
            processBuffer();
        }
        bafSitesSource.close();
        writer.close();
        return null;
    }

    /** Performs a Pearson's chi square test for goodness of fit to a biallelic model to determine
     * heterozygosity.  If this fails, BafEvidence.MISSING_VALUE is returned.  If the test succeeds,
     * the depth of the alt call as a fraction of the total depth is returned.
     */
    @VisibleForTesting double calcBAF( final LocusDepth locusDepth,
                                       final int refIndex, final int altIndex ) {
        final int totalDepth = locusDepth.getTotalDepth();
        if ( totalDepth < minTotalDepth ) {
            return BafEvidence.MISSING_VALUE;
        }
        final double expectRefAlt = totalDepth / 2.;
        final double altDepth = locusDepth.getDepth(altIndex);
        final double refDiff = locusDepth.getDepth(refIndex) - expectRefAlt;
        final double altDiff = altDepth - expectRefAlt;
        final double chiSq = (refDiff * refDiff + altDiff * altDiff) / expectRefAlt;
        final double fitProb = 1. - chiSqDist.cumulativeProbability(chiSq);
        if ( fitProb < minHetProbability ) {
            return BafEvidence.MISSING_VALUE;
        }
        return altDepth / totalDepth;
    }

    private boolean sameLocus( final Locatable locus ) {
        if ( sameLocusBuffer.isEmpty() ) {
            return true;
        }
        final LocusDepth curLocus = sameLocusBuffer.get(0);
        return curLocus.getContig().equals(locus.getContig()) &&
                curLocus.getStart() == locus.getStart();
    }

    private void processBuffer() {
        final VariantContext vc = nextSite();
        if ( !sameLocus(vc) ) {
            final Locatable loc = sameLocusBuffer.get(0);
            throw new UserException("expecting locus " + loc.getContig() + ":" + loc.getStart() +
                    ", but found locus " + vc.getContig() + ":" + vc.getStart() + " in baf sites vcf");
        }
        final List<Allele> alleles = vc.getAlleles();
        final int refIndex = Nucleotide.decode(alleles.get(0).getBases()[0]).ordinal();
        if ( refIndex > 3 ) {
            throw new UserException("ref call is not [ACGT] in vcf at " + vc.getContig() + ":" + vc.getStart());
        }
        final int altIndex = Nucleotide.decode(alleles.get(1).getBases()[0]).ordinal();
        if ( altIndex > 3 ) {
            throw new UserException("alt call is not [ACGT] in vcf at " + vc.getContig() + ":" + vc.getStart());
        }
        final List<BafEvidence> beList = new ArrayList<>(sameLocusBuffer.size());
        for ( final LocusDepth ld : sameLocusBuffer ) {
            final double baf = calcBAF(ld, refIndex, altIndex);
            if ( baf != BafEvidence.MISSING_VALUE ) {
                beList.add(new BafEvidence(ld.getSample(), ld.getContig(), ld.getStart(), baf));
            }
        }
        final int nBafs = beList.size();
        if ( nBafs == 1 ) {
            writer.write(new BafEvidence(beList.get(0), .5));
        } else {
            final RunningAverage ra = new RunningAverage();
            for ( final BafEvidence bafEvidence : beList ) {
                ra.add(bafEvidence.getValue());
            }
            final double stddev = ra.stddev();
            if ( stddev <= maxStdDev ) {
                final double[] vals = new double[nBafs];
                for ( int idx = 0; idx != nBafs; ++idx ) {
                    vals[idx] = beList.get(idx).getValue();
                }
                Arrays.sort(vals);
                final int midIdx = nBafs / 2;
                final double median =
                        (nBafs & 1) != 0 ? vals[midIdx] : (vals[midIdx] + vals[midIdx - 1])/2.;
                final double adj = .5 - median;
                for ( final BafEvidence bafEvidence : beList ) {
                    writer.write(new BafEvidence(bafEvidence, bafEvidence.getValue()+adj));
                }
            }
        }
        sameLocusBuffer.clear();
    }

    private VariantContext nextSite() {
        while ( bafSitesIterator.hasNext() ) {
            final VariantContext vc = bafSitesIterator.next();
            if ( vc.isSNP() && vc.isBiallelic() ) {
                return vc;
            }
        }
        throw new UserException("baf sites vcf exhausted before locus depth data");
    }
}
