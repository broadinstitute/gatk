package org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading;

import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.graphs.MultiSampleEdge;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public final class TestingReadThreadingGraph extends ReadThreadingGraph {
    /*************************************************************
     * Simple string representation support for testing purposes *
     *************************************************************/

    private static final Pattern PROPERTIES_PATTERN = Pattern.compile("^\\s*\\[[^\\]]*\\]");
    private static final Pattern PATH_PATTERN = Pattern.compile("\\{((\\S+):)?([^\\}]*)\\}");
    private static final Pattern KMERSIZE_EXTRACTOR_PATTERN = Pattern.compile("^\\s*\\[[^\\]]*(ks|kmerSize)\\s*=\\s*(\\d+)\\s*[,\\]]");
    private static final long serialVersionUID = 1l;


    /**
     * Constructs a read-threading-graph for a string representation.
     *
     * <p>
     *     Note: only used for testing.
     *     Checkout {@see HaplotypeGraphUnitTest} for examples.
     * </p>
     * @param s the string representation of the graph {@code null}.
     */
    public TestingReadThreadingGraph(final String s) {
        super(kmerSizeFromString(s));
        applyString(s);
        setAlreadyBuilt();
    }

    /**
     * Obtain the kmer size for the string representation.
     * @param str the source string representation.
     * @return 1 or greater.
     * @throws IllegalArgumentException if {@code} str does not contain a valid representation.
     */
    private static int kmerSizeFromString(final String str) {
        final Matcher matcher = KMERSIZE_EXTRACTOR_PATTERN.matcher(str);
        if (matcher.find()) {
            return Integer.parseInt(matcher.group(2));
        } else {
            throw new IllegalArgumentException("the input graph spec does not indicate the kmerSize");
        }
    }

    /**
     * Apply description string into the graph.
     *
     * <p>
     *     Note: this is done just for testing purposes.
     *     Checkout {@see HaplotypeGraphUnitTest} for examples.
     * </p>
     * @param str the string representation.
     */
    private void applyString(final String str) {
        final Matcher propertiesSectionMatcher = PROPERTIES_PATTERN.matcher(str);
        final int pathStart = propertiesSectionMatcher.find() ? propertiesSectionMatcher.end() : 0;

        final String pathString = str.substring(pathStart);
        final Matcher pathMatcher = PATH_PATTERN.matcher(pathString);

        boolean referenceFound = false;
        final Map<String,MultiDeBruijnVertex> vertexById = new HashMap<>();

        // Loop between path strings and add them one by one.
        while (pathMatcher.find()) {
            final String label = pathMatcher.group(2);
            final boolean isReference = "REF".equals(label);
            if (referenceFound) {
                if (isReference) {
                    throw new IllegalArgumentException("there are two reference paths");
                }
            } else if ( isReference ) {
                referenceFound = true;
            }

            // Divide each path into its elements getting a list of sequences and labels if applies:
            final String elementsString = pathMatcher.group(3);
            final String[] elements = elementsString.split("\\s*->\\s*");
            if (elements.length == 0) {
                throw new IllegalArgumentException("empty path not allowed");
            }
            final String[] seqs = new String[elements.length];
            final String[] ids = new String[elements.length];
            for (int i = 0; i < elements.length; i++) {
                ids[i] = pathElementId(elements[i]);
                seqs[i] = pathElementSeq(elements[i]);
                if (seqs[i].isEmpty() && ids[i] == null) {
                    throw new IllegalArgumentException("path with empty element without an id");
                }
            }
            final boolean isSource =  ids[0] == null || !vertexById.containsKey(ids[0]);
            if (isSource && seqs[0].length() != kmerSize) {
                throw new IllegalArgumentException("source sequence length must be the same as the kmerSize "
                        + ids[0] + ' ' + seqs[0] + ' ' + pathMatcher.group());
            }
            final MultiDeBruijnVertex firstVertex;
            if (ids[0] != null && vertexById.containsKey(ids[0])) {
                firstVertex = vertexById.get(ids[0]);
            } else {
                firstVertex = new MultiDeBruijnVertex(seqs[0].getBytes());
                addVertex(firstVertex);
                if (ids[0] != null) {
                    vertexById.put(ids[0], firstVertex);
                }
            }
            if (!seqs[0].isEmpty() &&
                    ((isSource && !firstVertex.getSequenceString().equals(seqs[0]))
                            || (!isSource && firstVertex.getSuffix() != seqs[0].getBytes()[0]))) {
                throw new IllegalArgumentException("mismatched first element sequence");
            }

            MultiDeBruijnVertex lastVertex = firstVertex;
            for (int i = 1; i < elements.length; i++) {
                if (seqs[i].length() > 1) {
                    throw new IllegalArgumentException("non-source vertex sequence must have length 1");
                }
                final MultiDeBruijnVertex nextVertex;
                if (ids[i] == null || !vertexById.containsKey(ids[i])) {
                    final Set<MultiDeBruijnVertex> nextVertices = getNextVertices(lastVertex,seqs[i].getBytes()[0]);
                    if (nextVertices.isEmpty()) {
                        nextVertex = new MultiDeBruijnVertex(extendSequence(lastVertex.getSequence(),seqs[i].getBytes()[0]));
                        addVertex(nextVertex);
                    } else {
                        nextVertex = nextVertices.iterator().next();
                    }
                    if (ids[i] != null) {
                        vertexById.put(ids[i], nextVertex);
                    }
                } else {
                    nextVertex = vertexById.get(ids[i]);
                }
                final MultiSampleEdge edge = addEdge(lastVertex,nextVertex);
                if (isReference) {
                    edge.setIsRef(true);
                }
                lastVertex = nextVertex;
            }
        }
    }

    private static String pathElementId(final String element) {
        final int openBracketPosition = element.indexOf('(');

        if (openBracketPosition == -1) {
            return null;
        }

        final int closeBracketPosition = element.lastIndexOf(')');
        if (closeBracketPosition == -1) {
            throw new IllegalArgumentException("non-closed id parantesys found in element: " + element);
        }
        final String result = element.substring(openBracketPosition + 1,closeBracketPosition).trim();
        if (result.isEmpty()) {
            throw new IllegalArgumentException("empty id found in element: " + element);
        }
        return result;
    }

    /**
     * Returns the lenght of a path element in the string representation.
     * @param element the query element.
     * @return 0 or greater.
     */
    private static String pathElementSeq(final String element) {
        final int parentesysPos = element.indexOf('(');

        if (parentesysPos == -1) {
            return element.trim();
        }

        return element.substring(0,parentesysPos).trim();
    }

    /**
     * Add a base to the end of a byte sequence.
     * @param sequence sequence where to add the base to.
     * @param b base to add.
     * @return never {@code null}, a new array each time.
     */
    private static byte[] extendSequence(final byte[] sequence, final byte b) {
        final byte[] result = new byte[sequence.length];
        System.arraycopy(sequence, 1, result, 0, sequence.length - 1);
        result[result.length - 1] = b;
        return result;
    }

    @Override
    public TestingReadThreadingGraph clone() {
        return (TestingReadThreadingGraph) super.clone();
    }
}
