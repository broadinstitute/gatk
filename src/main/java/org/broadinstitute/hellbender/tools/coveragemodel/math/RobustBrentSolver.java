package org.broadinstitute.hellbender.tools.coveragemodel.math;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.AbstractUnivariateSolver;
import org.apache.commons.math3.analysis.solvers.BrentSolver;
import org.apache.commons.math3.exception.NoBracketingException;
import org.apache.commons.math3.exception.TooManyEvaluationsException;
import org.apache.commons.math3.util.FastMath;
import org.broadinstitute.hellbender.utils.Utils;

import javax.annotation.Nullable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * A robust version of the Brent solver that tries to avoid spurious non-bracketing conditions.
 *
 * The root is bracketed by searching the solution search interval. If multiple roots are found and
 * a non-null merit function is provided, the root with maximum merit is chosen. If no merit function
 * is provided or if merits are equal, the leftmost root is chosen by convention.
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 */
public final class RobustBrentSolver extends AbstractUnivariateSolver {
    private final UnivariateFunction meritFunc;
    private final int numBisections;
    private final int depth;

    public RobustBrentSolver(final double relativeAccuracy,
                             final double absoluteAccuracy,
                             final double functionValueAccuracy,
                             @Nullable final UnivariateFunction meritFunc,
                             final int numBisections,
                             final int depth) {
        super(relativeAccuracy, absoluteAccuracy, functionValueAccuracy);
        this.meritFunc = meritFunc;
        this.numBisections = numBisections;
        this.depth = depth;
    }

    @Override
    public double solve(final int maxEval,
                        final UnivariateFunction objFunc,
                        final double min,
                        final double max)
            throws TooManyEvaluationsException, NoBracketingException {
        setup(maxEval, objFunc, min, max, min); /* the last parameter is not actually used */
        return doSolve();
    }

    @Override
    protected double doSolve() throws TooManyEvaluationsException, NoBracketingException {
        final double min = getMin();
        final double max = getMax();
        final double[] xSearchGrid = createHybridSearchGrid(min, max, numBisections, depth);
        final double[] fSearchGrid = Arrays.stream(xSearchGrid).map(this::computeObjectiveValue).toArray();

        /* find bracketing intervals on the search grid */
        final List<Bracket> bracketsList = detectBrackets(xSearchGrid, fSearchGrid);
        if (bracketsList.isEmpty()) {
            throw new NoBracketingException(min, max, fSearchGrid[0], fSearchGrid[fSearchGrid.length-1]);
        }
        final BrentSolver solver = new BrentSolver(getRelativeAccuracy(), getAbsoluteAccuracy(), getFunctionValueAccuracy());
        final List<Double> roots = bracketsList.stream()
                .map(b -> solver.solve(getMaxEvaluations(), this::computeObjectiveValue, b.min, b.max, 0.5 * (b.min + b.max)))
                .collect(Collectors.toList());
        if (roots.size() == 1 || meritFunc == null) {
            return roots.get(0);
        }
        final double[] merits = roots.stream().mapToDouble(meritFunc::value).toArray();
        final int bestRootIndex = IntStream.range(0, roots.size())
                .boxed()
                .max((i, j) -> (int) (merits[i] - merits[j]))
                .get();
        return roots.get(bestRootIndex);
    }

    /**
     * Generates a hybrid search grid concentrated around {@code min}. The hybrid grid starts with a base-2
     * logarithmic grid. Each grid element is further divided uniformly into {@code refinementOrder} intervals.
     *
     * @param min left endpoint
     * @param max right endpoint
     * @param logSubdivisions number of logarithmic refinements
     * @param uniformSubdivisions number of uniform refinements
     * @return a double array of grid points
     */
    @VisibleForTesting
    double[] createHybridSearchGrid(final double min, final double max, final int logSubdivisions, final int uniformSubdivisions) {
        final double[] baseGrid = createLogarithmicGrid(min, max, logSubdivisions, 2);
        final int nBaseGrid = logSubdivisions + 2;
        if (uniformSubdivisions > 1) {
            final double[] refinedGrid = new double[(nBaseGrid - 1) * uniformSubdivisions + 1];
            for (int j = 0; j < nBaseGrid - 1; j++) {
                final double len = (baseGrid[j + 1] - baseGrid[j]) / uniformSubdivisions;
                for (int k = 0; k < uniformSubdivisions; k++) {
                    refinedGrid[j * uniformSubdivisions + k] = baseGrid[j] + k * len;
                }
            }
            refinedGrid[(nBaseGrid - 1) * uniformSubdivisions] = max;
            return refinedGrid;
        } else {
            return baseGrid;
        }
    }

    /**
     * Creates a logarithmic grid concentrated around {@code min}.
     *
     * @param min left endpoint
     * @param max right endpoint
     * @param subdivisions number of logarithmic subdivisions
     * @param base logarithm base
     * @return a double array of grid points of length {@code subdivisions} + 2
     */
    private double[] createLogarithmicGrid(final double min, final double max, final int subdivisions, final double base) {
        Utils.validateArg(base > 1, "The logarithm base must be greater than 1");
        final double[] grid = new double[subdivisions + 2];
        grid[0] = 0;
        grid[subdivisions + 1] = max - min;
        for (int j = subdivisions; j > 0; j--) {
            grid[j] = grid[j+1] / base;
        }
        for (int j = 0; j < subdivisions + 2; j++) {
            grid[j] += min;
        }
        return grid;
    }

    @VisibleForTesting
    static List<Bracket> detectBrackets(final double[] x, final double[] f) {
        final List<Bracket> brackets = new ArrayList<>();
        final double[] signs = new double[f.length];
        for (int i = 0; i < f.length; i++) {
            signs[i] = FastMath.signum(f[i]);
        }
        double prevSignum = signs[0];
        int prevIdx = 0;
        int idx = 1;
        while (idx < f.length) {
            if (signs[idx]*prevSignum <= 0) {
                brackets.add(new Bracket(x[prevIdx], x[idx]));
                prevIdx = idx;
                prevSignum = signs[idx];
            }
            idx++;
        }
        return brackets;
    }

    @VisibleForTesting
    static final class Bracket {
        final double min, max;

        Bracket(final double min, final double max) {
            this.min = min;
            this.max = max;
        }
    }
}
