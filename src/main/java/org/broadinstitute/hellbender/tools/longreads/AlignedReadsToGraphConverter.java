package org.broadinstitute.hellbender.tools.longreads;

import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.esotericsoftware.kryo.serializers.JavaSerializer;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.argparser.ExperimentalFeature;
import org.broadinstitute.hellbender.cmdline.programgroups.CoverageAnalysisProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadWalker;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.longreads.graph.AlignedBaseGraphCollection;
import org.broadinstitute.hellbender.tools.longreads.graph.LabeledEdgeType;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.util.Collections;
import java.util.List;

@CommandLineProgramProperties(
        summary = "Make a graph from a bam file containing reads aligned to a reference.  This tool will always create a .gfa file representing the graph of the aligned reads.",
        oneLineSummary = "Make a graph from a bam file containing reads aligned to a reference.",
        programGroup = CoverageAnalysisProgramGroup.class
)
@ExperimentalFeature
public class AlignedReadsToGraphConverter extends ReadWalker {
    private static final Logger logger = LogManager.getLogger(AlignedReadsToGraphConverter.class);

    //==================================================================================================================
    // Public Static Members:

    public static final String OUTPUT_FILE_BASE_NAME_ARG = "output-file-base-name";
    public static final String SKIP_ZIP_ARG = "skip-zip";
    public static final String GRAPH_IN_ARG = "graph-in";
    public static final String GRAPH_OUT_ARG = "graph-out";
    public static final String RELABEL_EDGE_TYPES_ARG = "relabel-edge-types";
    public static final String CREATE_DOT_FILES_ARG = "create-dot-files";
    public static final String CREATE_GEXF_FILES_ARG = "create-gexf-files";
    public static final String EDGE_TYPE_ARG = "edge-type";

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Public Members:

    @Argument(
            fullName  = OUTPUT_FILE_BASE_NAME_ARG,
            optional = true,
            doc = "Base name for output files to be written.")
    private String outputFileBaseName = "aligned_base_graph";

    @Advanced
    @Argument(
            fullName  = SKIP_ZIP_ARG,
            optional = true,
            doc = "Skip the zipping stage in graph creation.")
    private Boolean skipZip = false;

    @Advanced
    @Argument(
            fullName  = GRAPH_IN_ARG,
            optional = true,
            mutex = {"edge-type"},
            doc = "Initialize the graph object with the file at given path.")
    private File inputGraphFile = null;

    @Advanced
    @Argument(
            fullName  = GRAPH_OUT_ARG,
            optional = true,
            doc = "Saves the graph to the given path.")
    private File outputGraphFile = null;

    @Advanced
    @Argument(
            fullName  = RELABEL_EDGE_TYPES_ARG,
            optional = true,
            doc = "If true, will overwrite types of existing edges when new data are added.")
    private boolean relabelEdgeTypes = false;

    @Argument(
            fullName  = CREATE_DOT_FILES_ARG,
            optional = true,
            doc = "Create an additional DOT file for each GFA file created.")
    private boolean createDotFiles = false;

    @Argument(
            fullName  = CREATE_GEXF_FILES_ARG,
            optional = true,
            doc = "Create an additional GEXF file for each GFA file created.\n" +
                  "WARNING: GEXF files can get VERY large.  Use this at your own peril.")
    private Boolean createGexfFiles = false;

    @Advanced
    @Argument(
            fullName  = EDGE_TYPE_ARG,
            optional = true,
            mutex = {"graph-in"},
            doc = "Use a specific edge type when creating links in the graph.  Will only work when creating a graph from scratch.")
    private LabeledEdgeType edgeType = LabeledEdgeType.LABELED_EDGE;

    //==================================================================================================================
    // Private Members:

    private AlignedBaseGraphCollection alignedBaseGraphCollection;

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Collections.emptyList();
    }

    @Override
    public void onTraversalStart() {
        if ( inputGraphFile == null ) {
            alignedBaseGraphCollection = new AlignedBaseGraphCollection(edgeType);
        }
        else {
            alignedBaseGraphCollection = deserializeGraphFromFile(inputGraphFile);
        }

        // Propagate our edge relabeling:
        if ( relabelEdgeTypes ) {
            logger.info("Relabeling edge types is ENABLED.");
        }
        alignedBaseGraphCollection.setRelabelEdgeTypes(relabelEdgeTypes);
    }

    @Override
    public void apply(final GATKRead read, final ReferenceContext referenceContext, final FeatureContext featureContext) {
        // Add in our sequence data here:
        alignedBaseGraphCollection.addSequence( read );
    }

    @Override
    public Object onTraversalSuccess() {
        if ( !skipZip ) {
            // Make sure we zip all adjacent nodes into big boy nodes because they're all grown up now.
            alignedBaseGraphCollection.collapseAdjacentNodes();
        }

        if ( createDotFiles ) {
            logger.info("Writing DOT files...");
            alignedBaseGraphCollection.serializeToDotFiles(outputFileBaseName);
        }

        if ( createGexfFiles ) {
            logger.info("Writing GEXF files...");
            alignedBaseGraphCollection.serializeToGexfFiles(outputFileBaseName);
        }

        // Write our master GFA files:
        logger.info("Writing GFA 1 files...");
        alignedBaseGraphCollection.serializeToGfa1Files(outputFileBaseName);

        // DISABLED until we determine if this is actually valid.
//        logger.info("Writing GFA 2 files...");
//        alignedBaseGraphCollection.serializeToGfa2Files(outputFileBaseName);

        // Save the graph if we should save it:
        if ( outputGraphFile != null ) {
            serializeGraphToFile(alignedBaseGraphCollection, outputGraphFile);
        }

        return null;
    }

    //==================================================================================================================
    // Static Methods:

//    /**
//     * Serialize graph data to a file to use in later runs.
//     * @param alignedBaseGraphCollection {@link AlignedBaseGraphCollection} containing graph data to serialize to a file.
//     * @param outputGraphFile {@link File} to which to write serialized graph data.
//     */
//    private static void serializeGraphToFile(final AlignedBaseGraphCollection alignedBaseGraphCollection,
//                                             final File outputGraphFile) {
//        logger.info("Writing graph to file (" + outputGraphFile.getAbsolutePath() + ") ...");
//
//        try( final FileOutputStream fileOutputStream = new FileOutputStream(outputGraphFile);
//            final ObjectOutputStream out = new ObjectOutputStream(fileOutputStream)) {
//            out.writeObject(alignedBaseGraphCollection);
//        }
//        catch ( final FileNotFoundException ex ) {
//            throw new UserException("Can't save graph to output file.", ex);
//        }
//        catch ( final IOException ex ) {
//            throw new GATKException("Can't save graph to output file.", ex);
//        }
//    }
//
//    /**
//     * Deserialize graph data from a file to initialize our graph collection.
//     * @param inputGraphFile {@link File} from which to read serialized graph data.
//     * @return A {@link AlignedBaseGraphCollection} containing the data represented in the given {@code inputGraphFile}.
//     */
//    private static AlignedBaseGraphCollection deserializeGraphFromFile(final File inputGraphFile) {
//        logger.info("Initializing graph from file (" + inputGraphFile.getAbsolutePath() + ") ...");
//
//        final AlignedBaseGraphCollection alignedBaseGraphCollection;
//        try (final FileInputStream fileInputStream = new FileInputStream(inputGraphFile);
//             final ObjectInputStream inputStream = new ObjectInputStream(fileInputStream)) {
//
//            alignedBaseGraphCollection = (AlignedBaseGraphCollection) inputStream.readObject();
//        }
//        catch ( final FileNotFoundException ex ) {
//            throw new UserException("Can't open graph input file.", ex);
//        }
//        catch ( final Exception ex ) {
//            throw new GATKException("WTF!?!?!?", ex);
//        }
//
//        return alignedBaseGraphCollection;
//    }

    /**
     * Serialize graph data to a file to use in later runs.
     * @param alignedBaseGraphCollection {@link AlignedBaseGraphCollection} containing graph data to serialize to a file.
     * @param outputGraphFile {@link File} to which to write serialized graph data.
     */
    private static void serializeGraphToFile(final AlignedBaseGraphCollection alignedBaseGraphCollection,
                                             final File outputGraphFile) {
        logger.info("Writing graph to kryo file (" + outputGraphFile.getAbsolutePath() + ")...");
        final Kryo kryo = new Kryo();
        kryo.register(AlignedBaseGraphCollection.class, new JavaSerializer());

        try (final Output output = new Output(new FileOutputStream(outputGraphFile))) {
            kryo.writeObject(output, alignedBaseGraphCollection);
        }
        catch ( final FileNotFoundException ex ) {
            throw new UserException("Can't save graph to output file.", ex);
        }
    }

    /**
     * Deserialize graph data from a file to initialize our graph collection.
     * @param inputGraphFile {@link File} from which to read serialized graph data.
     * @return A {@link AlignedBaseGraphCollection} containing the data represented in the given {@code inputGraphFile}.
     */
    private static AlignedBaseGraphCollection deserializeGraphFromFile(final File inputGraphFile) {
        logger.info("Initializing graph from kryo file (" + inputGraphFile.getAbsolutePath() + ") ...");
        final Kryo kryo = new Kryo();
        kryo.register(AlignedBaseGraphCollection.class, new JavaSerializer());

        final AlignedBaseGraphCollection alignedBaseGraphCollection;
        try (final Input input = new Input(new FileInputStream(inputGraphFile))) {
            alignedBaseGraphCollection = kryo.readObject(input, AlignedBaseGraphCollection.class);
        }
        catch ( final FileNotFoundException ex ) {
            throw new UserException("Can't open graph input file.", ex);
        }
        catch ( final Exception ex ) {
            throw new GATKException("WTF!?!?!?", ex);
        }

        logger.info("Read in graph collection with reads: " + alignedBaseGraphCollection.getNumSequencesAdded());

        return alignedBaseGraphCollection;
    }

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
