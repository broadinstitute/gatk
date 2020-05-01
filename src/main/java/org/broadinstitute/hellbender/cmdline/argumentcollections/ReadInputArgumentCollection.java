package org.broadinstitute.hellbender.cmdline.argumentcollections;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.GATKPathSpecifier;
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
                    "for each input will be inferred automatically.",
            common = true,
            optional = true)
    protected List<GATKPathSpecifier> readIndices;

    //Lazily computed the first time it is requested
    private List<ReadIndexPair> readIndexPairs = null;

    public static class ReadIndexPair{
        private final GATKPathSpecifier reads;
        private final GATKPathSpecifier index;

        public ReadIndexPair(final GATKPathSpecifier reads, final GATKPathSpecifier index) {
            this.reads = reads;
            this.index = index;
        }

        public GATKPathSpecifier getReads() {
            return reads;
        }

        public GATKPathSpecifier getIndex() {
            return index;
        }
    }

    /**
     * Get the raw list of BAM/SAM/CRAM files specified at the command line.
     * Paths are the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     * These will be processed to resolve bundle files.
     */
    protected abstract List<GATKPathSpecifier> getRawReadPathSpecifiers();


    public List<GATKPathSpecifier> getReadPathSpecifiers(){
        return getReadIndexPairs().stream().map(ReadIndexPair::getReads).collect(Collectors.toList());
    }

    public List<ReadIndexPair> getReadIndexPairs() {
        //check if it's already been cached
        if( readIndexPairs == null){
            //compute it if necessary
            final List<GATKPathSpecifier> rawReadPathSpecifiers = getRawReadPathSpecifiers();
            final int numberOfReadSourcesSpecified = rawReadPathSpecifiers.size();
            final int numberOfReadIndexesSpecified = readIndices.size();
            if( !readIndices.isEmpty() && numberOfReadSourcesSpecified != numberOfReadIndexesSpecified){
                throw new UserException("If  --"+ StandardArgumentDefinitions.READ_INDEX_LONG_NAME + " is specified " +
                        "it must be specified once for every read input that is specified. " +
                        "\n Found " + numberOfReadSourcesSpecified +"  read sources and " + numberOfReadIndexesSpecified + " read indexes.");
            }
            final List<ReadIndexPair> pairs = new ArrayList<>(numberOfReadIndexesSpecified);
            for( int i = 0; i < numberOfReadSourcesSpecified ; i++){
                //This has the problem where we can't identify a .json that doesn't have the right extension
                final GATKPathSpecifier rawReadPath = rawReadPathSpecifiers.get(i);
                final ReadIndexPair pair;
                if(ReadsBundle.looksLikeReadsBundle(rawReadPath)){
                    if( !readIndices.isEmpty()){
                        throw new UserException("You can specify read/index pairs with json read bundles " +
                                "OR with the --"+ StandardArgumentDefinitions.READ_INDEX_LONG_NAME+ " argument but you cannot mix the two.");
                    }
                    final ReadsBundle readsBundle = ReadsBundle.fromPath(rawReadPath);
                    pair = new ReadIndexPair(readsBundle.getReads(), readsBundle.getIndex());
                } else {
                    pair = new ReadIndexPair(rawReadPath, readIndices.get(i));
                }
                pairs.add(pair);
            }
            readIndexPairs = Collections.unmodifiableList(pairs);
        }
        return readIndexPairs;
    }

    /**
     * Get the list of BAM/SAM/CRAM files specified at the command line.
     * Paths are the preferred format, as this can handle both local disk and NIO direct access to cloud storage.
     */
    public List<Path> getReadPaths() {
        return getReadPathSpecifiers().stream().map(GATKPathSpecifier::toPath).collect(Collectors.toList());
    }

    /**
     * @return The list of indices to be used with the read inputs, or {@code null} if none were specified and the indices should be
     *         inferred automatically.
     *
     *         If explicit indices are specified, they must be specified for all read inputs, and are assumed to be in the same
     *         order as the read inputs.
     */
    public List<Path> getReadIndexPaths() {
        if ( readIndices == null || readIndices.isEmpty() ) {
            return null;
        }
        return getReadIndexPairs().stream()
                .map(ReadIndexPair::getIndex)
                .map(GATKPathSpecifier::toPath)
                .collect(Collectors.toList());
    }

    /**
     * Get the read validation stringency specified at the command line, or the default value if none was specified
     * at the command line.
     */
    public ValidationStringency getReadValidationStringency() { return readValidationStringency; };
}
