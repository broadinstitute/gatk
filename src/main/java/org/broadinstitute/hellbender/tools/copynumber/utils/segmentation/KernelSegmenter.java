package org.broadinstitute.hellbender.tools.copynumber.utils.segmentation;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.DefaultRealMatrixChangingVisitor;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.SingularValueDecomposition;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.RandomGeneratorFactory;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.utils.optimization.PersistenceOptimizer;
import org.broadinstitute.hellbender.utils.IndexRange;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * <p>
 * Segments data (i.e., finds multiple changepoints) using a method based on the kernel-segmentation algorithm
 * described in <a href="https://hal.inria.fr/hal-01413230/document">https://hal.inria.fr/hal-01413230/document</a>,
 * which gives a framework to quickly calculate the cost of a segment given a low-rank approximation to a specified kernel.
 * However, unlike the algorithm described there, which seeks to minimize a global segmentation cost,
 * the method implemented here instead finds candidate changepoints based on the local costs within
 * windows of various sizes.
 * </p>
 *
 * <p>
 * Given <i>N</i> sequential data points of type {@code DATA} to segment, the basic steps of the method are:
 *
 * <ol>
 *     1) Select <i>C<sub>max</sub></i>, the maximum number of changepoints to discover.
 * </ol>
 * <ol>
 *     2) Select a kernel (linear for sensitivity to changes in the distribution mean,
 *     Gaussian with a specified variance for multimodal data, etc.)
 *     and a subsample of <i>p</i> points to approximate it using singular value decomposition.
 *     This is provided via a lambda of the form {@code (DATA, DATA) -> double}.
 * </ol>
 * <ol>
 *     3) Select window sizes <i>w<sub>j</sub></i> for which to calculate local costs at each point.
 *     To be precise, we calculate the cost of a changepoint at the point with index <i>i</i>,
 *     assuming adjacent segments containing the points with indices <i>[i - w<sub>j</sub> + 1, i]</i>
 *     and <i>[i + 1, i + w<sub>j</sub>]</i>.
 * </ol>
 * <ol>
 *     4) For each of these cost functions, find (up to) the  <i>C<sub>max</sub></i> most significant local minima.
 *     The problem of finding local minima of a noisy function can be solved by using topological persistence
 *     (e.g., <a href="https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html">https://people.mpi-inf.mpg.de/~weinkauf/notes/persistence1d.html</a>
 *     and <a href="http://www2.iap.fr/users/sousbie/web/html/indexd3dd.html?post/Persistence-and-simplification">http://www2.iap.fr/users/sousbie/web/html/indexd3dd.html?post/Persistence-and-simplification</a>).
 * </ol>
 * <ol>
 *     5) These sets of local minima from all window sizes together provide the pool of candidate changepoints
 *     (some of which may overlap exactly or approximately).  We perform backward selection using the global segmentation cost.
 *     That is, we calculate the global segmentation cost given all the candidate changepoints,
 *     calculate the cost change for removing each of the changepoints individually,
 *     remove the changepoint with the minimum cost change, and repeat.
 *     This gives the global cost as a function of the number of changepoints <i>C</i>.
 * </ol>
 * <ol>
 *     6) Add a penalty <i>A * C + B * C * log (N / C)</i> to the global cost and find the minimum to determine the
 *     number of changepoints, where <i>A</i> and <i>B</i> are specified penalty factors.
 *
 * </ol>
 * </p>
 *
 * <p>
 * See discussion at <a href="https://github.com/broadinstitute/gatk/issues/2858#issuecomment-324125586">https://github.com/broadinstitute/gatk/issues/2858#issuecomment-324125586</a>
 * and accompanying plots for more detail.
 * </p>
 *
 * <p>
 * Note that we break with camelCase naming convention in places to match some notation in the paper
 * </p>
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class KernelSegmenter<DATA> {
    public enum ChangepointSortOrder {
        BACKWARD_SELECTION, INDEX
    }

    private static final Logger logger = LogManager.getLogger(KernelSegmenter.class);

    private static final int RANDOM_SEED = 1216;
    private static final double EPSILON = 1E-10;

    private final List<DATA> data;

    public KernelSegmenter(final List<DATA> data) {
        this.data = Collections.unmodifiableList(new ArrayList<>(Utils.nonNull(data)));
    }

    /**
     * Returns a list of the indices of the changepoints, either sorted by decreasing change to the global segmentation cost
     * or by increasing index order.
     * @param maxNumChangepoints                    maximum number of changepoints to return (first and last points do not count towards this number)
     * @param kernel                                kernel function used to calculate segment costs
     * @param kernelApproximationDimension          dimension of low-rank approximation to the kernel
     * @param windowSizes                           list of sizes to use for the flanking segments used to calculate local changepoint costs
     * @param numChangepointsPenaltyLinearFactor    factor A for penalty of the form A * C, where C is the number of changepoints
     * @param numChangepointsPenaltyLogLinearFactor factor B for penalty of the form B * C * log (N / C),
     *                                              where C is the number of changepoints and N is the number of data points
     * @param changepointSortOrder                  sort by decreasing change to the global segmentation cost or by increasing index order
     */
    public List<Integer> findChangepoints(final int maxNumChangepoints,
                                          final BiFunction<DATA, DATA, Double> kernel,
                                          final int kernelApproximationDimension,
                                          final List<Integer> windowSizes,
                                          final double numChangepointsPenaltyLinearFactor,
                                          final double numChangepointsPenaltyLogLinearFactor,
                                          final ChangepointSortOrder changepointSortOrder) {
        ParamUtils.isPositiveOrZero(maxNumChangepoints, "Maximum number of changepoints must be non-negative.");
        ParamUtils.isPositive(kernelApproximationDimension, "Dimension of kernel approximation must be positive.");
        Utils.validateArg(!windowSizes.isEmpty(), "At least one window size must be provided.");
        Utils.validateArg(windowSizes.stream().allMatch(ws -> ws > 0), "Window sizes must all be positive.");
        Utils.validateArg(windowSizes.stream().distinct().count() == windowSizes.size(), "Window sizes must all be unique.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLinearFactor,
                "Linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");
        ParamUtils.isPositiveOrZero(numChangepointsPenaltyLogLinearFactor,
                "Log-linear factor for the penalty on the number of changepoints per chromosome must be non-negative.");

        if (maxNumChangepoints == 0) {
            logger.warn("No changepoints were requested, returning an empty list...");
            return Collections.emptyList();
        }
        if (data.isEmpty()) {
            logger.warn("No data points were provided, returning an empty list...");
            return Collections.emptyList();
        }

        logger.info(String.format("Finding up to %d changepoints in %d data points...", maxNumChangepoints, data.size()));
        final RandomGenerator rng = RandomGeneratorFactory.createRandomGenerator(new Random(RANDOM_SEED));

        logger.info("Calculating low-rank approximation to kernel matrix...");
        final RealMatrix reducedObservationMatrix = calculateReducedObservationMatrix(rng, data, kernel, kernelApproximationDimension);
        final double[] kernelApproximationDiagonal = calculateKernelApproximationDiagonal(reducedObservationMatrix);

        logger.info(String.format("Finding changepoint candidates for all window sizes %s...", windowSizes.toString()));
        final List<Integer> changepointCandidates = findChangepointCandidates(
                data, reducedObservationMatrix, kernelApproximationDiagonal, maxNumChangepoints, windowSizes);

        logger.info("Performing backward model selection on changepoint candidates...");
        return selectChangepoints(
                changepointCandidates, maxNumChangepoints, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor,
                reducedObservationMatrix, kernelApproximationDiagonal).stream()
                .sorted((a, b) -> changepointSortOrder.equals(ChangepointSortOrder.INDEX) ? Integer.compare(a, b) : 0)    //if BACKWARD_SELECTION, simply retain original order from backward model selection
                .collect(Collectors.toList());
    }

    private static final class Segment {
        private final int start;    //inclusive index of start point
        private final int end;      //inclusive index of end point
        private final double cost;

        private Segment(final int start,
                        final int end,
                        final double cost) {
            this.start = start;
            this.end = end;
            this.cost = cost;
        }

        private Segment(final int start,
                        final int end,
                        final RealMatrix reducedObservationMatrix,
                        final double[] kernelApproximationDiagonal) {
            this(start, end, calculateSegmentCost(start, end, reducedObservationMatrix, kernelApproximationDiagonal).C);
        }
    }

    //represents some quantities used to calculate segment costs iteratively
    private static final class Cost {
        private final double D;     //diagonal term
        private final double[] W;   //intermediate quantity for calculating off-diagonal terms
        private final double V;     //intermediate quantity for calculating off-diagonal terms
        private final double C;     //total cost

        private Cost(final double D,
                     final double[] W,
                     final double V,
                     final double C) {
            this.D = D;
            this.W = W;
            this.V = V;
            this.C = C;
        }
    }

    //calculates the N x p reduced observation matrix, defined as Z in equation preceding Eq. 14 in https://hal.inria.fr/hal-01413230/document
    private static <DATA> RealMatrix calculateReducedObservationMatrix(final RandomGenerator rng,
                                                                       final List<DATA> data,
                                                                       final BiFunction<DATA, DATA, Double> kernel,
                                                                       final int kernelApproximationDimension) {
        if (kernelApproximationDimension > data.size()) {
            logger.warn("Specified dimension of the kernel approximation exceeds the number of data points to segment; " +
                    "using all data points to calculate kernel matrix.");
        }

        //subsample data with replacement
        final int numSubsample = Math.min(kernelApproximationDimension, data.size());
        logger.info(String.format("Subsampling %d points from data to find kernel approximation...", numSubsample));
        final List<DATA> dataSubsample = numSubsample == data.size()
                ? data
                : IntStream.range(0, numSubsample).mapToObj(i -> data.get(rng.nextInt(data.size()))).collect(Collectors.toList());

        //calculate (symmetric) kernel matrix of subsampled data
        logger.info(String.format("Calculating kernel matrix of subsampled data (%d x %d)...", numSubsample, numSubsample));
        final RealMatrix subKernelMatrix = new Array2DRowRealMatrix(numSubsample, numSubsample);
        for (int i = 0; i < numSubsample; i++) {
            for (int j = 0; j < i; j++) {
                final double value = kernel.apply(dataSubsample.get(i), dataSubsample.get(j));
                subKernelMatrix.setEntry(i, j, value);
                subKernelMatrix.setEntry(j, i, value);
            }
            subKernelMatrix.setEntry(i, i, kernel.apply(dataSubsample.get(i), dataSubsample.get(i)));
        }

        //perform SVD of kernel matrix of subsampled data
        logger.info(String.format("Performing SVD of kernel matrix of subsampled data (%d x %d)...", numSubsample, numSubsample));
        final SingularValueDecomposition svd = new SingularValueDecomposition(subKernelMatrix);

        //calculate reduced observation matrix
        logger.info(String.format("Calculating reduced observation matrix (%d x %d)...", data.size(), numSubsample));
        final double[] invSqrtSingularValues = Arrays.stream(svd.getSingularValues()).map(Math::sqrt).map(x -> 1. / (x + EPSILON)).toArray();
        final RealMatrix subKernelUMatrix = new Array2DRowRealMatrix(numSubsample, numSubsample);
        subKernelUMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int i, int j, double value) {
                return svd.getU().getEntry(i, j) * invSqrtSingularValues[j];
            }
        });
        final RealMatrix reducedKernelMatrix = new Array2DRowRealMatrix(data.size(), numSubsample);
        reducedKernelMatrix.walkInOptimizedOrder(new DefaultRealMatrixChangingVisitor() {
            @Override
            public double visit(int i, int j, double value) {
                return kernel.apply(data.get(i), dataSubsample.get(j));
            }
        });
        return reducedKernelMatrix.multiply(subKernelUMatrix);
    }

    //for N x p matrix Z_ij, returns the N-dimensional vector sum(Z_ij * Z_ij, j = 0,..., p - 1),
    //which are the diagonal elements K_ii of the approximate kernel matrix
    private static double[] calculateKernelApproximationDiagonal(final RealMatrix reducedObservationMatrix) {
        return new IndexRange(0, reducedObservationMatrix.getRowDimension())
                .mapToDouble(i -> MathUtils.square(reducedObservationMatrix.getRowVector(i).getNorm()));
    }

    //finds indices of changepoint candidates from all window sizes
    private static <DATA> List<Integer> findChangepointCandidates(final List<DATA> data,
                                                                  final RealMatrix reducedObservationMatrix,
                                                                  final double[] kernelApproximationDiagonal,
                                                                  final int maxNumChangepoints,
                                                                  final List<Integer> windowSizes) {
        final List<Integer> changepointCandidates = new ArrayList<>(windowSizes.size() * maxNumChangepoints);

        //for each window size, calculate local changepoint costs at each point and add maxNumChangepoints candidates
        //(this is overkill, but we cannot guarantee that the most significant maxNumChangepoints changepoints
        //do not all appear at only a single window size)
        for (final int windowSize : windowSizes) {
            logger.debug(String.format("Calculating local changepoints costs for window size %d...", windowSize));
            if (windowSize > data.size()) {
                logger.warn(String.format("Number of points needed to calculate local changepoint costs (2 * window size = %d) " +
                        "exceeds number of data points %d.  Local changepoint costs will not be calculated for this window size.",
                        2 * windowSize, data.size()));
                continue;
            }
            final double[] windowCosts = calculateWindowCosts(reducedObservationMatrix, kernelApproximationDiagonal, windowSize);

            logger.debug(String.format("Finding local minima of local changepoint costs for window size %d...", windowSize));
            final List<Integer> windowCostLocalMinima = new ArrayList<>(new PersistenceOptimizer(windowCosts).getMinimaIndices());
            windowCostLocalMinima.remove(Integer.valueOf(0));                //remove first data point if present
            windowCostLocalMinima.remove(Integer.valueOf(data.size() - 1));  //remove last data point if present
            changepointCandidates.addAll(windowCostLocalMinima.subList(0, Math.min(maxNumChangepoints, windowCostLocalMinima.size())));
        }

        if (changepointCandidates.isEmpty()) {
            throw new GATKException.ShouldNeverReachHereException("No changepoint candidates found.");
        }

        return changepointCandidates;
    }

    //performs backward model selection to order changepoints by increasing change to the global segmentation cost
    //and returns the requested number
    private static List<Integer> selectChangepoints(final List<Integer> changepointCandidates,
                                                    final int maxNumChangepoints,
                                                    final double numChangepointsPenaltyLinearFactor,
                                                    final double numChangepointsPenaltyLogLinearFactor,
                                                    final RealMatrix reducedObservationMatrix,
                                                    final double[] kernelApproximationDiagonal) {
        final List<Integer> changepoints = new ArrayList<>(changepointCandidates.size());

        //calculate penalties as a function of the number of changepoints
        final int numData = reducedObservationMatrix.getRowDimension();
        final List<Double> changepointPenalties = IntStream.range(0, maxNumChangepoints + 1)
                .mapToObj(numChangepoints -> calculateChangepointPenalty(
                        numChangepoints, numChangepointsPenaltyLinearFactor, numChangepointsPenaltyLogLinearFactor, numData))
                .collect(Collectors.toList());

        //construct initial list of all segments and initialize costs
        final List<Integer> candidateStarts = changepointCandidates.stream().sorted().distinct()
                .map(i -> Math.min(i + 1, numData - 1)).collect(Collectors.toList());
        candidateStarts.add(0, 0);
        final List<Integer> candidateEnds = changepointCandidates.stream().sorted().distinct().collect(Collectors.toList());
        candidateEnds.add(numData - 1);
        final int numSegments = candidateStarts.size();
        final List<Segment> segments = IntStream.range(0, numSegments)
                .mapToObj(i -> new Segment(candidateStarts.get(i), candidateEnds.get(i), reducedObservationMatrix, kernelApproximationDiagonal))
                .collect(Collectors.toList());
        final List<Double> totalSegmentationCosts = new ArrayList<>(Collections.singletonList(segments.stream().mapToDouble(s -> s.cost).sum()));
        final List<Double> costsForSegmentPairs = IntStream.range(0, numSegments - 1)
                .mapToObj(i -> segments.get(i).cost + segments.get(i + 1).cost)
                .collect(Collectors.toList());  //sum of the costs for the segments in each adjacent pair
        final List<Double> costsForMergedSegmentPairs = IntStream.range(0, numSegments - 1)
                .mapToObj(i -> new Segment(candidateStarts.get(i), candidateEnds.get(i + 1), reducedObservationMatrix, kernelApproximationDiagonal).cost)
                .collect(Collectors.toList());  //cost of each adjacent pair when considered as a single segment
        final List<Double> costsForMergingSegmentPairs = IntStream.range(0, numSegments - 1)
                .mapToObj(i -> costsForSegmentPairs.get(i) - costsForMergedSegmentPairs.get(i))
                .collect(Collectors.toList());  //cost for merging each adjacent pair into a single segment

        //iteratively merge the segment pair with greatest merge cost and update all costs until only a single segment remains
        for (int i = 0; i < numSegments - 1; i++) {
            //find segment pair to merge and calculate quantities for resulting merged segment
            final int indexOfLeftSegmentToMerge = costsForMergingSegmentPairs.indexOf(Collections.max(costsForMergingSegmentPairs));
            final double newCost = costsForMergedSegmentPairs.get(indexOfLeftSegmentToMerge);
            final int newStart = segments.get(indexOfLeftSegmentToMerge).start;
            final int mergepoint = segments.get(indexOfLeftSegmentToMerge).end;
            final int newEnd = segments.get(indexOfLeftSegmentToMerge + 1).end;

            //remove segment pair and insert merged segment into list of segments
            segments.remove(indexOfLeftSegmentToMerge);
            segments.remove(indexOfLeftSegmentToMerge);
            segments.add(indexOfLeftSegmentToMerge, new Segment(newStart, newEnd, newCost));

            //update segment-pair quantities
            costsForSegmentPairs.remove(indexOfLeftSegmentToMerge);
            costsForMergedSegmentPairs.remove(indexOfLeftSegmentToMerge);
            costsForMergingSegmentPairs.remove(indexOfLeftSegmentToMerge);
            if (indexOfLeftSegmentToMerge > 0) {                    //if segment pair that was merged was not the first pair, update segment-pair quantities using segment to left
                costsForSegmentPairs.set(indexOfLeftSegmentToMerge - 1, segments.get(indexOfLeftSegmentToMerge - 1).cost + segments.get(indexOfLeftSegmentToMerge).cost);
                costsForMergedSegmentPairs.set(indexOfLeftSegmentToMerge - 1, new Segment(segments.get(indexOfLeftSegmentToMerge - 1).start, newEnd, reducedObservationMatrix, kernelApproximationDiagonal).cost);
                costsForMergingSegmentPairs.set(indexOfLeftSegmentToMerge - 1, costsForSegmentPairs.get(indexOfLeftSegmentToMerge - 1) - costsForMergedSegmentPairs.get(indexOfLeftSegmentToMerge - 1));
            }
            if (indexOfLeftSegmentToMerge < segments.size() - 1) {  //if segment pair that was merged was not the last pair, update segment-pair quantities using segment to right
                costsForSegmentPairs.set(indexOfLeftSegmentToMerge, segments.get(indexOfLeftSegmentToMerge).cost + segments.get(indexOfLeftSegmentToMerge + 1).cost);
                costsForMergedSegmentPairs.set(indexOfLeftSegmentToMerge, new Segment(newStart, segments.get(indexOfLeftSegmentToMerge + 1).end, reducedObservationMatrix, kernelApproximationDiagonal).cost);
                costsForMergingSegmentPairs.set(indexOfLeftSegmentToMerge, costsForSegmentPairs.get(indexOfLeftSegmentToMerge) - costsForMergedSegmentPairs.get(indexOfLeftSegmentToMerge));
            }

            //update total segmentation costs and changepoints
            totalSegmentationCosts.add(0, segments.stream().mapToDouble(s -> s.cost).sum());
            changepoints.add(0, mergepoint);
        }

        //find optimal number of changepoints according to penalty function
        final int effectiveMaxNumChangepoints = Math.min(maxNumChangepoints, changepoints.size());
        final List<Double> totalSegmentationCostsPlusPenalties = IntStream.range(0, effectiveMaxNumChangepoints + 1)
                .mapToObj(i -> totalSegmentationCosts.get(i) + changepointPenalties.get(i))
                .collect(Collectors.toList());
        final int numChangepointsOptimal = totalSegmentationCostsPlusPenalties.indexOf(Collections.min(totalSegmentationCostsPlusPenalties));

        logger.info(String.format("Found %d changepoints after applying the changepoint penalty.", numChangepointsOptimal));
        return changepoints.subList(0, numChangepointsOptimal);
    }

    private static double calculateChangepointPenalty(final int numChangepoints,
                                                      final double numChangepointsPenaltyLinearFactor,
                                                      final double numChangepointsPenaltyLogLinearFactor,
                                                      final int numData) {
        return numChangepointsPenaltyLinearFactor * numChangepoints
                + numChangepointsPenaltyLogLinearFactor * numChangepoints * Math.log(numData / (numChangepoints + EPSILON));
    }

    /**
     * Calculates the cost of a segment.  This is defined by Eq. 11 of
     * <a href="https://hal.inria.fr/hal-01413230/document">https://hal.inria.fr/hal-01413230/document</a>
     * (except we use the low-rank approximation to the kernel, as described in Sec. 3.2, ibid).
     * Various recurrence relations are used to calculate costs iteratively.
     * @param start inclusive start index of segment
     * @param end   inclusive end index of segment
     * @param reducedObservationMatrix      N x p matrix of projected observations, where N is the number of data points
     *                                      and p is the dimension of the low-rank approximation of the kernel matrix;
     *                                      this is the Z matrix described in the text preceding Eq. 14, ibid
     * @param kernelApproximationDiagonal   N diagonal terms of the low-rank approximation to the kernel matrix
     */
    private static Cost calculateSegmentCost(final int start,
                                             final int end,
                                             final RealMatrix reducedObservationMatrix,
                                             final double[] kernelApproximationDiagonal) {
        final int N = reducedObservationMatrix.getRowDimension();
        final int p = reducedObservationMatrix.getColumnDimension();

        //initialize quantities for recurrence
        double D = kernelApproximationDiagonal[start];
        final double[] W = Arrays.copyOf(reducedObservationMatrix.getRow(start), p);
        double V = Arrays.stream(W).map(w -> w * w).sum();

        //generate indices for iteration; we need to wrap around to beginning of data if start > end
        final List<Integer> indices = start <= end 
                ? IntStream.range(start + 1, end + 1).boxed().collect(Collectors.toList()) 
                : IntStream.concat(IntStream.range(start + 1, N), IntStream.range(0, end + 1)).boxed().collect(Collectors.toList());

        //use recurrence relations to iteratively calculate cost
        for (final int tauPrime : indices) {
            D += kernelApproximationDiagonal[tauPrime];
            double ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(tauPrime, j) * W[j];
                W[j] += reducedObservationMatrix.getEntry(tauPrime, j);
            }
            V += 2. * ZdotW + kernelApproximationDiagonal[tauPrime];
        }
        final double C = D - V / (indices.size() + 1);

        return new Cost(D, W, V, C);
    }

    /**
     * Calculates the local costs at each point for a given window size <i>w</i>.  Using Eq. 11 of
     * <a href="https://hal.inria.fr/hal-01413230/document">https://hal.inria.fr/hal-01413230/document</a>
     * (except we use the low-rank approximation to the kernel, as described in Sec. 3.2, ibid), for each point
     * indexed by <i>i</i>, we calculate the cost of it being a changepoint with two flanking segments that
     * contain the points with indices <i>[i - w + 1, i]</i> and <i>[i + 1, i + w]</i>, respectively, and
     * subtract the cost of a single segment containing all of these points.
     * Various recurrence relations are used to calculate costs iteratively.
     * @param reducedObservationMatrix      N x p matrix of projected observations, where N is the number of data points
     *                                      and p is the dimension of the low-rank approximation of the kernel matrix;
     *                                      this is the Z matrix described in the text preceding Eq. 14, ibid
     * @param kernelApproximationDiagonal   N diagonal terms of the low-rank approximation to the kernel matrix
     * @param windowSize                    number of points to include in either flanking segment when calculating cost
     */
    private static double[] calculateWindowCosts(final RealMatrix reducedObservationMatrix,
                                                 final double[] kernelApproximationDiagonal,
                                                 final int windowSize) {
        final int N = reducedObservationMatrix.getRowDimension();
        final int p = reducedObservationMatrix.getColumnDimension();

        //initialize indices of the boundaries of the two flanking segments, wrapping around to beginning of data if necessary
        int center = 0;
        int start = (center - windowSize + 1 + N) % N;
        int end = (center + windowSize) % N;

        //initialize costs of flanking segments and total segment
        final Cost leftCost = calculateSegmentCost(start, center, reducedObservationMatrix, kernelApproximationDiagonal);
        final Cost rightCost = calculateSegmentCost(center + 1, end, reducedObservationMatrix, kernelApproximationDiagonal);
        final Cost totalCost = calculateSegmentCost(start, end, reducedObservationMatrix, kernelApproximationDiagonal);

        //initialize quantities for recurrence
        double leftD = leftCost.D;
        final double[] leftW = Arrays.copyOf(leftCost.W, p);
        double leftV = leftCost.V;
        double leftC = leftCost.C;

        double rightD = rightCost.D;
        final double[] rightW = Arrays.copyOf(rightCost.W, p);
        double rightV = rightCost.V;
        double rightC = rightCost.C;

        double totalD = totalCost.D;
        final double[] totalW = Arrays.copyOf(totalCost.W, p);
        double totalV = totalCost.V;
        double totalC = totalCost.C;

        final double[] windowCosts = new double[N];
        windowCosts[center] = leftC + rightC - totalC;

        double ZdotW;
        final double windowSizeReciprocal = 1. / windowSize;

        //slide segments along data and use recurrence relations to iteratively update costs
        for (center = 0; center < N; center++) {
            final int centerNext = (center + 1) % N;
            final int endNext = (end + 1) % N;

            //update quantities in left segment
            leftD -= kernelApproximationDiagonal[start];
            ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(start, j) * leftW[j];
                leftW[j] -= reducedObservationMatrix.getEntry(start, j);
            }
            leftV += -2. * ZdotW + kernelApproximationDiagonal[start];

            leftD += kernelApproximationDiagonal[centerNext];
            ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(centerNext, j) * leftW[j];
                leftW[j] += reducedObservationMatrix.getEntry(centerNext, j);
            }
            leftV += 2. * ZdotW + kernelApproximationDiagonal[centerNext];

            leftC = leftD - leftV * windowSizeReciprocal;

            //update quantities in right segment
            rightD -= kernelApproximationDiagonal[centerNext];
            ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(centerNext, j) * rightW[j];
                rightW[j] -= reducedObservationMatrix.getEntry(centerNext, j);
            }
            rightV += -2. * ZdotW + kernelApproximationDiagonal[centerNext];

            rightD += kernelApproximationDiagonal[endNext];
            ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(endNext, j) * rightW[j];
                rightW[j] += reducedObservationMatrix.getEntry(endNext, j);
            }
            rightV += 2. * ZdotW + kernelApproximationDiagonal[endNext];

            rightC = rightD - rightV * windowSizeReciprocal;

            //update quantities in total segment
            totalD -= kernelApproximationDiagonal[start];
            ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(start, j) * totalW[j];
                totalW[j] -= reducedObservationMatrix.getEntry(start, j);
            }
            totalV += -2. * ZdotW + kernelApproximationDiagonal[start];

            totalD += kernelApproximationDiagonal[endNext];
            ZdotW = 0.;
            for (int j = 0; j < p; j++) {
                ZdotW += reducedObservationMatrix.getEntry(endNext, j) * totalW[j];
                totalW[j] += reducedObservationMatrix.getEntry(endNext, j);
            }
            totalV += 2. * ZdotW + kernelApproximationDiagonal[endNext];

            totalC = totalD - 0.5 * totalV * windowSizeReciprocal;

            //record cost of changepoint at this position
            windowCosts[centerNext] = leftC + rightC - totalC;

            //slide windows
            start = (start + 1) % N;
            end = endNext;
        }
        return windowCosts;
    }
}