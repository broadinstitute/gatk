package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jgrapht.Graph;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.SimpleGraph;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Class to work with exclusive pairs of elements, example - pairs of alleles that do not occur in the haplotypes
 * @param <R> rows of matrix (e.g. haplotypes)
 * @param <C> columns of matrix (e.g. alleles)
 *
 * The class keeps a binary matrix of rows x columns and looks for columns that never have 1 in the same row. (e.g.
 *           alleles that never occur in the same haplotype.
 *           It then allows to generate a graph with edges: a non-coocurence relationship and find connected components
 *           in ths graph. This is useful for AlleleFiltering where we are basing on clusters of exclusive alleles.
 *
 */
public class OccurrenceMatrix<R,C> {
    protected static final Logger logger = LogManager.getLogger(OccurrenceMatrix.class);

    private List<R> rowNames;
    private List<C> colNames;
    private Map<C, Integer> col2idx;
    private int nRows;
    private int nCols;
    private boolean[][] occurrenceMatrix;

    /**
     *
     * @param input_map input appearance map between R -> C (haplotype to alleles that it contains)
     */
    public OccurrenceMatrix(final Map<R, Collection<C>> input_map){
        nRows = input_map.size();
        rowNames = input_map.keySet().stream().collect(Collectors.toList());

        col2idx = new HashMap<>();
        int col_count = 0;
        for ( R row: input_map.keySet() ){
            for (C col: input_map.get(row)){
                if (!col2idx.containsKey(col)){
                    col2idx.put(col, col_count);
                    col_count++;
                }
            }
        }

        nCols = col_count;
        colNames = new ArrayList<>();
        for (int i = 0 ; i < nCols; i++ ){
            colNames.add(null);
        }
        for ( C col: col2idx.keySet() ){
            colNames.set(col2idx.get(col), col);
        }


        occurrenceMatrix = new boolean[nRows][nCols];
        for (int r = 0; r < nRows; r++){
            for ( C col: input_map.get(rowNames.get(r)))
                occurrenceMatrix[r][col2idx.get(col)] = true;
        }
    }

    /**
     * find pairs of columns which do not both have a true value at any row
     * @return - pairs of found columns
     */
    public List<Pair<C, C>> nonCoOcurringColumns() {

        List<Pair<Integer, Integer>> result = new ArrayList<>();
        for (int i = 0; i < nCols; i++) {
            for (int j = i + 1; j < nCols; j++) {
                boolean flag = false;
                for (int r = 0; r < nRows; r++) {
                    if (occurrenceMatrix[r][i] & occurrenceMatrix[r][j]) {
                        flag = true;
                        break;
                    }
                }
                if (!flag) {
                    result.add(ImmutablePair.of(i,j));
                }

            }
        }
        List<Pair<C,C>> vc_result = new ArrayList<>();
        for (Pair<Integer, Integer> res: result) {
            vc_result.add(ImmutablePair.of(colNames.get(res.getLeft()), colNames.get(res.getRight())));
        }
        return vc_result;
    }

    /**
     * Analyse columns for being independent of each other.
     * @param nonCoOcurringColumns - pairs of columns which do not both have a true value at any row (non-cooccurring)
     * @return - list of columns sets. The columns in each set are non-cooccurring with at least one other member
     */
    public List<Set<C>> getIndependentSets(List<Pair<C,C>> nonCoOcurringColumns){
        Graph<C, DefaultEdge> nonConnectedAllelesGraph = new SimpleGraph<>(DefaultEdge.class);
        colNames.stream().forEach(x -> nonConnectedAllelesGraph.addVertex(x));
        nonCoOcurringColumns.stream().forEach(edge->nonConnectedAllelesGraph.addEdge(edge.getLeft(), edge.getRight()));

        ConnectivityInspector<C, DefaultEdge> ci = new ConnectivityInspector<>(nonConnectedAllelesGraph);
        List<Set<C>> result = ci.connectedSets();

        // debug log messages
        if ( logger.isDebugEnabled() ) {
            logger.debug(String.format("GIS: Received %d alleles that generate %d connected components", colNames.size(), result.size()));
            logger.debug("GIS: Here are the components:");

            for (int i = 0; i < result.size(); i++) {
                String str = new String();
                for (C allele : result.get(i)) {
                    str += " ";
                    str += allele.toString();
                }
                logger.debug(String.format("---- GIS: (%d) %s", i, str));
            }
        }

        return result;
    }

}