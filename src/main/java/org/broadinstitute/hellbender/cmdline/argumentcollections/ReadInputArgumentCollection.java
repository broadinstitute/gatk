package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPath;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.bundle.ReadsBundle;
import org.broadinstitute.hellbender.utils.read.ReadConstants;

import java.io.Serializable;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;


/**
 * An abstract argument collection for use with tools that accept input files containing reads
 * (eg., BAM/SAM/CRAM files).
 */
public abstract class ReadInputArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_LONG_NAME,
            shortName = StandardArgumentDefinitions.READ_VALIDATION_STRINGENCY_SHORT_NAME,
            doc = "Validation stringency for all SAM/BAM/CRAM/SRA files read by this program.  The default stringency value SILENT " +
                    "can improve performance when processing a BAM file in which variable-length data (read, qualities, tags) " +
                    "do not otherwise need to be decoded.",
            common=true,
            optional=true)
    protected ValidationStringency readValidationStringency = ReadConstants.DEFAULT_READ_VALIDATION_STRINGENCY;


    @Argument(fullName = StandardArgumentDefinitions.READ_INDEX_LONG_NAME, shortName = StandardArgumentDefinitions.READ_INDEX_SHORT_NAME,
            doc = "Indices to use for the read inputs. If specified, an index must be provided for every read input " +
                    "and in the same order as the read inputs. If this argument is not specified, the path to the index " +
                    "for each input will be inferred automatically.  This argument is deprecated.  Use a read-bundle.json instead.",
            common = true,
            optional = true)
    @Deprecated
    protected List<GATKPath> readIndices;

    @Advanced
    @Argument(fullName = "dont-infer-bam-indexes", doc = "Don't attempt to infer the location of bam indexes based on the filename.")
    protected boolean dontInferBamIndexes = false;

    //Lazily computed the first time it is requested
    private List<ReadIndexPair> readIndexPairs = null;

    /**
     * Get the raw list of BAM/SAM/CRAM inputs specified at the command line.
     * GATKPath is the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     * These will be processed to resolve bundle files.
     */
    protected abstract List<GATKPath> getRawReadPathSpecifiers();


    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line.
     * Paths are the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */

    public List<Path> getReadPaths() {
        return getReadPathSpecifiers().stream().map(GATKPath::toPath).collect(Collectors.toList());
    }

    /**
     *
     * @return
     */
    public List<GATKPath> getReadPathSpecifiers(){
        return getReadIndexPairs().stream().map(ReadIndexPair::getReads).collect(Collectors.toList());
    }


    /**
     * Get the matched pairs of BAM/SAM/CRAM and resolved indexes that were specified on the command line
     * @return
     */
    public List<ReadIndexPair> getReadIndexPairs() {
        //check if it's already been cached
        if( readIndexPairs == null){
            //compute it if necessary
            final List<GATKPath> rawReadPathSpecifiers = getRawReadPathSpecifiers();
            final int numberOfReadSourcesSpecified = rawReadPathSpecifiers.size();
            final int numberOfReadIndexesSpecified = readIndices.size();
            if( !readIndices.isEmpty() && numberOfReadSourcesSpecified != numberOfReadIndexesSpecified){
                throw new UserException("If  --"+ StandardArgumentDefinitions.READ_INDEX_LONG_NAME + " is specified " +
                        "it must be specified once for every read input that is specified. " +
                        "\n Found " + numberOfReadSourcesSpecified +"  read sources and " + numberOfReadIndexesSpecified + " read indexes.");
            }

            if ( !readIndices.isEmpty() && rawReadPathSpecifiers.stream().anyMatch(ReadsBundle::looksLikeAReadsBundle)){
                if( !readIndices.isEmpty()){
                    throw new UserException("You can specify read/index pairs with json read bundles " +
                            "OR with the --"+ StandardArgumentDefinitions.READ_INDEX_LONG_NAME+ " argument but you cannot mix the two.");
                }
            }

            final List<ReadIndexPair> pairs = new ArrayList<>(numberOfReadIndexesSpecified);
            for( int i = 0; i < numberOfReadSourcesSpecified ; i++){
                //TODO This has the problem where we can't identify a .json that doesn't have the right extension
                final GATKPath rawReadPath = rawReadPathSpecifiers.get(i);
                final ReadIndexPair pair;
                if(ReadsBundle.looksLikeAReadsBundle(rawReadPath)){
                    //if it looks like a bundle, load it
                    final ReadsBundle readsBundle = ReadsBundle.fromPath(rawReadPath);
                    pair = new ReadIndexPair(readsBundle.getReads(), readsBundle.getIndex());
                } else if (!readIndices.isEmpty()) {
                    //if it isn't a bundle and we have read indexes provided than get it from the list
                    pair = new ReadIndexPair(rawReadPath, readIndices.get(i));
                } else {
                    //otherwise we have to decide to infer the index path or not
                    if(dontInferBamIndexes) {
                        //in this case we explicitly set the index to null since it wasn't specified
                        pair = new ReadIndexPair(rawReadPath, null);
                    } else {
                        //otherwise we try to guess
                        pair = ReadIndexPair.guessPairFromReads(rawReadPath);
                    }
                }
                pairs.add(pair);
            }
            readIndexPairs = Collections.unmodifiableList(pairs);
        }
        return readIndexPairs;
    }

    /**
     * Get the read validation stringency specified at the command line, or the default value if none was specified
     * at the command line.
     */
    public ValidationStringency getReadValidationStringency() { return readValidationStringency; };

}
