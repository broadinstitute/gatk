package org.broadinstitute.hellbender.tools.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.tribble.Feature;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.engine.MultiFeatureWalker;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.MathUtils.RunningAverage;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodec;
import org.broadinstitute.hellbender.utils.codecs.FeatureOutputCodecFinder;
import org.broadinstitute.hellbender.utils.codecs.FeatureSink;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * <p>Merges locus-sorted LocusDepth evidence files, and calculates the bi-allelic frequency for each
 * sample and site, and writes these values as a BafEvidence output file.</p>
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
    public static final String SAMPLE_NAMES_NAME = "sample-names";
    public static final String BAF_EVIDENCE_FILE_NAME = "baf-evidence-output";
    public static final String MIN_TOTAL_DEPTH_NAME = "min-total-depth";
    public static final String MAX_STD_NAME = "max-std";
    public static final String MIN_HET_LIKELIHOOD = "min-het-likelihood";
    public static final String COMPRESSION_LEVEL_NAME = "compression-level";
    public static ChiSquaredDistribution chiSqDist = new ChiSquaredDistribution(null, 1.);

    @Argument(
            doc = "Locus depth files to process",
            fullName = LOCUS_DEPTH_FILE_NAME,
            shortName = StandardArgumentDefinitions.FEATURE_SHORT_NAME )
    private List<GATKPath> locusDepthFiles;

    @Argument(fullName = SAMPLE_NAMES_NAME, doc = "Sample names", optional = true)
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
            doc = "Minimum likelihood of site being biallelic and heterozygous",
            fullName = MIN_HET_LIKELIHOOD, optional = true )
    private double minHetLikelihood = .5;

    @Argument(
            doc = "Output compression level",
            fullName = COMPRESSION_LEVEL_NAME,
            minValue = 0, maxValue = 9, optional = true )
    private int compressionLevel = 4;

    private FeatureSink<BafEvidence> writer;
    private final List<BafEvidence> sameLocusBuffer = new ArrayList<>();

    @Override
    @SuppressWarnings("unchecked")
    public void onTraversalStart() {
        super.onTraversalStart();
        final FeatureOutputCodec<? extends Feature, ? extends FeatureSink<? extends Feature>> codec =
                FeatureOutputCodecFinder.find(bafOutputFile);
        if ( !codec.getFeatureType().isAssignableFrom(BafEvidence.class) ) {
            throw new UserException("We're intending to write BafEvidence, but the feature type " +
                    "associated with the output file expects features of type " +
                    codec.getFeatureType().getSimpleName());
        }
        writer = (FeatureSink<BafEvidence>)codec.makeSink(bafOutputFile,
                                    getBestAvailableSequenceDictionary(),
                                    sampleNames, compressionLevel);
    }

    @Override
    public void apply( LocusDepth locusDepth, Object header,
                       ReadsContext readsContext, ReferenceContext referenceContext ) {
        final double baf = calcBAF(locusDepth);
        if ( baf > 0. ) {
            final BafEvidence bafEvidence =
                    new BafEvidence(locusDepth.getSample(), locusDepth.getContig(),
                                    locusDepth.getStart(), baf);
            if ( !sameLocus(bafEvidence) ) {
                processBuffer();
            }
            sameLocusBuffer.add(bafEvidence);
        }
    }

    @Override
    public Object onTraversalSuccess() {
        super.onTraversalSuccess();
        writer.close();
        return null;
    }

    @VisibleForTesting double calcBAF( final LocusDepth locusDepth ) {
        final int totalDepth = locusDepth.getTotalDepth();
        if ( totalDepth < minTotalDepth ) {
            return 0.;
        }
        double expectRefAlt = totalDepth / 2.;
        final double refDiff = locusDepth.getRefDepth() - expectRefAlt;
        final double altDiff = locusDepth.getAltDepth() - expectRefAlt;
        double chiSq = (refDiff * refDiff + altDiff * altDiff) / expectRefAlt;
        double fitProb = 1. - chiSqDist.cumulativeProbability(chiSq);
        if ( fitProb < minHetLikelihood ) {
            return 0.;
        }
        return (double)locusDepth.getAltDepth() / totalDepth;
    }

    private boolean sameLocus( final BafEvidence bafEvidence ) {
        if ( sameLocusBuffer.isEmpty() ) {
            return true;
        }
        final BafEvidence curLocus = sameLocusBuffer.get(0);
        return curLocus.getContig().equals(bafEvidence.getContig()) &&
                curLocus.getStart() == bafEvidence.getStart();
    }

    private void processBuffer() {
        int nBafs = sameLocusBuffer.size();
        if ( nBafs == 1 ) {
            writer.write(new BafEvidence(sameLocusBuffer.get(0), .5));
        } else {
            final RunningAverage ra = new RunningAverage();
            for ( final BafEvidence bafEvidence : sameLocusBuffer ) {
                ra.add(bafEvidence.getValue());
            }
            final double stddev = ra.stddev();
            if ( stddev <= maxStdDev ) {
                final double[] vals = new double[nBafs];
                for ( int idx = 0; idx != nBafs; ++idx ) {
                    vals[idx] = sameLocusBuffer.get(idx).getValue();
                }
                Arrays.sort(vals);
                int midIdx = nBafs / 2;
                double median =
                        (nBafs & 1) != 0 ? vals[midIdx] : (vals[midIdx] + vals[midIdx - 1])/2.;
                final double adj = .5 - median;
                for ( final BafEvidence bafEvidence : sameLocusBuffer ) {
                    writer.write(new BafEvidence(bafEvidence, bafEvidence.getValue()+adj));
                }
            }
        }
        sameLocusBuffer.clear();
    }
}
