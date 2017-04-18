package org.broadinstitute.hellbender.tools.coveragemodel.math;

import org.apache.commons.math3.analysis.solvers.AbstractUnivariateSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * This class implements a synchronized univariate solver for solving multiple independent equations.
 * It is to be used in situations where function queries have a costly overhead, though, simultaneous
 * queries of multiple functions have the same overhead as single queries.
 *
 * Consider the tasking of solving N independent equations:
 *
 *      f_1(x_1) = 0,
 *      f_2(x_2) = 0,
 *      ...
 *      f_N(x_N) = 0
 *
 * One approach is to solve these equations sequentially. In certain situations, each function
 * evaluation may be cheap but could entail a costly overhead (e.g. if the functions are evaluated
 * in a distributed architecture). It is desirable to minimize this overhead by bundling as many
 * function calls as possible, and querying the function in "chunks".
 *
 * Consider the ideal situation where function evaluations are infinitely cheap, however, each query has
 * a considerable overhead time of \tau. Also, let us assume that the overhead of simultaneously
 * querying {f_1(x_1), ..., f_N(x_N)} is the same as that of a single query, i.e. f_i(x_i). If the
 * univariate solver requires k queries on average, the overhead cost of the sequential approach is
 * O(k N \tau). By making simultaneous queries, this class reduces the overhead to O(k \tau).
 *
 * This is achieved by instantiating N threads for the N univariate solvers, accumulating their queries
 * and suspending them until all threads announce their required query.
 *
 * TODO github/gatk-protected issue #853
 * @implNote In the current implementation, we make a thread for each solver. This is fine if the number
 * of equations is reasonably small (< 200). In the future, the class must take a max number
 * of threads and limit concurrency.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class SynchronizedUnivariateSolver {

    /**
     * Default value for the absolute accuracy of function evaluations
     */
    private static final double DEFAULT_FUNCTION_ACCURACY = 1e-15;

    /**
     * Stores queries from instantiated solvers
     */
    private final ConcurrentHashMap<Integer, Double> queries;

    /**
     * Stores function evaluations on {@link #queries}
     */
    private final ConcurrentHashMap<Integer, Double> results;

    /**
     * Number of queries before making a function call
     */
    private final int numberOfQueriesBeforeCalling;

    /**
     * The objective functions
     */
    private final Function<Map<Integer, Double>, Map<Integer, Double>> func;

    /**
     * A list of solver jobs
     */
    private final List<UnivariateSolverJobDescription> jobDescriptions;
    private final List<UnivariateSolverSpecifications> solverDescriptions;
    private final Function<UnivariateSolverSpecifications, AbstractUnivariateSolver> solverFactory;
    private final Set<Integer> jobIndices;

    private final Lock resultsLock = new ReentrantLock();
    private final Condition resultsAvailable = resultsLock.newCondition();
    private CountDownLatch solversCountDownLatch;

    /**
     * Public constructor for invoking one of {@link AbstractUnivariateSolver}
     *
     * @param func the objective function (must be able to evaluate multiple univariate functions in one query)
     * @param numberOfQueriesBeforeCalling Number of queries before making a function call (the default value is
     *                                     the number of equations)
     */
    public SynchronizedUnivariateSolver(final Function<Map<Integer, Double>, Map<Integer, Double>> func,
                                        final Function<UnivariateSolverSpecifications, AbstractUnivariateSolver> solverFactory,
                                        final int numberOfQueriesBeforeCalling) {
        this.func = Utils.nonNull(func);
        this.solverFactory = Utils.nonNull(solverFactory);
        this.numberOfQueriesBeforeCalling = ParamUtils.isPositive(numberOfQueriesBeforeCalling, "Number of queries" +
                " before calling function evaluations must be positive");

        queries = new ConcurrentHashMap<>(numberOfQueriesBeforeCalling);
        results = new ConcurrentHashMap<>(numberOfQueriesBeforeCalling);
        jobDescriptions = new ArrayList<>();
        solverDescriptions = new ArrayList<>();
        jobIndices = new HashSet<>();
    }

    /**
     * Add a solver jobDescription
     *
     * @param index a unique index for the equation
     * @param min lower bound of the root
     * @param max upper bound of the root
     * @param x0 initial guess
     * @param absoluteAccuracy absolute accuracy
     * @param relativeAccuracy relative accuracy
     * @param functionValueAccuracy function value accuracy
     * @param maxEval maximum number of allowed evaluations
     */
    public void add(final int index, final double min, final double max, final double x0,
               final double absoluteAccuracy, final double relativeAccuracy,
               final double functionValueAccuracy, final int maxEval) {
        if (jobIndices.contains(index)) {
            throw new IllegalArgumentException("A jobDescription with index " + index + " already exists; jobDescription indices must" +
                    " be unique");
        }
        if (x0 <= min || x0 >= max) {
            throw new IllegalArgumentException(String.format("The initial guess \"%f\" for equation number \"%d\" is" +
                    " must lie inside the provided search bracket [%f, %f]", x0, index, min, max));
        }
        jobDescriptions.add(new UnivariateSolverJobDescription(index, min, max, x0, maxEval));
        solverDescriptions.add(new UnivariateSolverSpecifications(absoluteAccuracy, relativeAccuracy, functionValueAccuracy));
    }

    /**
     * Add a solver jobDescription using the default function accuracy {@link #DEFAULT_FUNCTION_ACCURACY}
     *
     * @param index a unique index for the equation
     * @param min lower bound of the root
     * @param max upper bound of the root
     * @param x0 initial guess
     * @param absoluteAccuracy absolute accuracy
     * @param relativeAccuracy relative accuracy
     * @param maxEval maximum number of allowed evaluations
     */
    public void add(final int index, final double min, final double max, final double x0,
                    final double absoluteAccuracy, final double relativeAccuracy,
                    final int maxEval) {
        add(index, min, max, x0, absoluteAccuracy, relativeAccuracy, DEFAULT_FUNCTION_ACCURACY, maxEval);
    }

    /**
     * Solve the equations
     *
     * @return a map from equation indices to the summary of results
     * @throws InterruptedException if any of the solver threads are interrupted
     */
    public Map<Integer, UnivariateSolverSummary> solve() throws InterruptedException {
        if (jobDescriptions.isEmpty()) {
            return Collections.emptyMap();
        }
        final Map<Integer, SolverWorker> solvers = new HashMap<>(jobDescriptions.size());
        solversCountDownLatch = new CountDownLatch(jobDescriptions.size());
        IntStream.range(0, jobDescriptions.size())
                .forEach(jobIdx -> solvers.put(jobDescriptions.get(jobIdx).getIndex(),
                        new SolverWorker(solverDescriptions.get(jobIdx), jobDescriptions.get(jobIdx))));

        /* start solver threads */
        solvers.values().forEach(worker -> new Thread(worker).start());

        /* wait for all workers to finish */
        solversCountDownLatch.await();

        return solvers.entrySet()
                .stream()
                .collect(Collectors.toMap(Map.Entry::getKey, entry -> entry.getValue().getSummary()));
    }

    /**
     * Require an evaluation of equation {@param index} at {$param x}
     *
     * @param index equation index
     * @param x equation argument
     * @return evaluated function value
     * @throws InterruptedException if the waiting thread is interrupted
     */
    private double evaluate(final int index, final double x) throws InterruptedException {
        queries.put(index, x);
        resultsLock.lock();
        final double value;
        try {
            fetchResults();
            while (!results.containsKey(index)) {
                resultsAvailable.await();
            }
            value = results.get(index);
            results.remove(index);
        } finally {
            resultsLock.unlock();
        }
        return value;
    }

    /**
     * Check if enough queries are in. If so, make a call and signal the waiting threads
     */
    private void fetchResults() {
        resultsLock.lock();
        try {
            if (queries.size() >= FastMath.min(numberOfQueriesBeforeCalling, solversCountDownLatch.getCount())) {
                results.putAll(func.apply(queries));
                queries.clear();
                resultsAvailable.signalAll();
            }
        } finally {
            resultsLock.unlock();
        }
    }

    public enum UnivariateSolverStatus {
        /**
         * Solution could not be bracketed
         */
        NO_BRACKETING,

        /**
         * Too many function evaluations
         */
        TOO_MANY_EVALUATIONS,

        /**
         * The solver found the root successfully
         */
        SUCCESS,

        /**
         * The status is not determined yet
         */
        TBD
    }

    /**
     * Stores the summary of a univariate solver jobDescription
     */
    public final class UnivariateSolverSummary {
        public final double x;
        public final int evaluations;
        public final UnivariateSolverStatus status;

        UnivariateSolverSummary(final double x, final int evaluations, final UnivariateSolverStatus status) {
            this.x = x;
            this.evaluations = evaluations;
            this.status = status;
        }
    }

    /**
     * A runnable version a solver
     */
    private final class SolverWorker implements Runnable {
        final UnivariateSolverJobDescription jobDescription;
        final UnivariateSolverSpecifications solverDescription;
        UnivariateSolverStatus status;
        private UnivariateSolverSummary summary;

        SolverWorker(final UnivariateSolverSpecifications solverDescription,
                     final UnivariateSolverJobDescription jobDescription) {
            this.solverDescription = solverDescription;
            this.jobDescription = jobDescription;
            status = UnivariateSolverStatus.TBD;
        }

        @Override
        public void run() {
            double sol;
            final AbstractUnivariateSolver abstractSolver = solverFactory.apply(solverDescription);
            try {
                sol = abstractSolver.solve(jobDescription.getMaxEvaluations(), x -> {
                    final double value;
                    try {
                        value = evaluate(jobDescription.getIndex(), x);
                    } catch (final InterruptedException ex) {
                        throw new RuntimeException(String.format("Evaluation of equation (n=%d) was interrupted --" +
                                " can not continue", jobDescription.getIndex()));
                    }
                    return value;
                }, jobDescription.getMin(), jobDescription.getMax(), jobDescription.getInitialGuess());
                status = UnivariateSolverStatus.SUCCESS;
            } catch (final NoBracketingException ex) {
                status = UnivariateSolverStatus.NO_BRACKETING;
                sol = Double.NaN;
            } catch (final TooManyEvaluationsException ex) {
                status = UnivariateSolverStatus.TOO_MANY_EVALUATIONS;
                sol = Double.NaN;
            }
                    summary = new UnivariateSolverSummary(sol, abstractSolver.getEvaluations(), status);
            solversCountDownLatch.countDown();
            fetchResults();
        }

        UnivariateSolverSummary getSummary() {
            return Utils.nonNull(summary, "Solver summary is not available");
        }
    }
}
