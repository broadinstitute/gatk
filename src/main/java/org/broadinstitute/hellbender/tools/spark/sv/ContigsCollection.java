package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import org.apache.commons.lang3.StringUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.exceptions.GATKException;
import scala.Tuple2;

import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Represents a collection of assembled contigs (not including the variants) produced by local assembler,
 * e.g. sga or fermi-lite.
 */
@VisibleForTesting
final class ContigsCollection implements Serializable {
    private static final long serialVersionUID = 1L;
    private static final Logger logger = LogManager.getLogger(ContigsCollection.class);

    private final List<LocalAssemblyContig> assemblyContents;

    public List<LocalAssemblyContig> getContents(){
        return assemblyContents;
    }

    @Override
    public String toString(){
        return StringUtils.join(toListOfStrings(),"\n");
    }

    /**
     * Pack the entire fasta record on one line so that it's easy for downstream
     * Spark tools to parse without worrying about partitioning. If the ContigCollection
     * is empty this returns null.
     */
    String toPackedFasta(){
        return StringUtils.join(toListOfStrings(),"|");
    }

    /**
     * To ease file format difficulties
     * @param packedFastaLine
     */
    @VisibleForTesting
    protected static ContigsCollection fromPackedFasta(final long assemblyID, final String packedFastaLine) {
        final List<String> fileContents = Arrays.asList(packedFastaLine.split("\\|"));
        if (fileContents.size() % 2 != 0) {
            throw new GATKException("Odd number of lines in breakpoint fasta" + packedFastaLine);
        }
        return new ContigsCollection(assemblyID, fileContents);
    }

    /**
     * @return {@code null} when {@link #assemblyContents} is,
     *          otherwise a list of string representation of contig id, and contig sequence
     */
    List<String> toListOfStrings(){
        if(null==assemblyContents){
            return null;
        }
        return assemblyContents.stream().flatMap(contig -> Stream.of(contig.contigID, contig.seq)).collect(Collectors.toList());
    }

    /**
     * Ctor to be used for SGA output
     */
    ContigsCollection(final long assemblyID, final List<String> sgaFASTAFileContents){

        if(null==sgaFASTAFileContents){
            assemblyContents = null;
        }else{
            assemblyContents = new ArrayList<>();
            for(int i=0; i<sgaFASTAFileContents.size(); i+=2){
                assemblyContents.add( new LocalAssemblyContig(assemblyID, sgaFASTAFileContents.get(i), sgaFASTAFileContents.get(i+1)) );
            }
        }
    }

    /**
     * Loads an RDD of {@link ContigsCollection} objects keyed by assembly ID from disk. The input file
     * should be the output of as {@link RunSGAViaProcessBuilderOnSpark}.
     */
    static JavaPairRDD<String, ContigsCollection> loadContigsCollectionKeyedByAssemblyId(final JavaSparkContext ctx, final String inputPath) {
        final JavaRDD<String> inputAssemblies = ctx.textFile(inputPath).cache();

        final JavaPairRDD<String, String> contigCollectionByAssemblyId =
                inputAssemblies
                        .flatMapToPair(ContigsCollection::splitAssemblyLine);

        return contigCollectionByAssemblyId.mapToPair(p -> new Tuple2<>(p._1(), fromPackedFasta(Long.parseLong(p._1().replace(SVConstants.FASTQ_OUT_PREFIX, "")), p._2())));
    }

    /**
     * input format is the text representation of an alignment region
     * @param alignedAssembledContigLine An input line with the tab-separated fields of an alignment region
     * @return A tuple with the breakpoint ID and string representation of an BreakpointAlignment, or an empty iterator if the line did not have two comma-separated values
     */
    static AlignmentRegion parseAlignedAssembledContigLine(final String alignedAssembledContigLine) {
        final String[] split = alignedAssembledContigLine.split("\t", -1);
        return AlignmentRegion.fromString(split);
    }

    /**
     * input format is tab separated AssemblyID, PackedFastaFile
     * @param assemblyLine An input line with a breakpoint ID and packed FASTA line
     * @return A tuple with the assembly ID and packed FASTA line, or an empty iterator if the line did not have two tab-separated values
     */
    static Iterable<Tuple2<String, String>> splitAssemblyLine(final String assemblyLine) {

        final String[] split = assemblyLine.split("\t");
        if (split.length < 2) {
            logger.info("No assembled contigs for breakpoint " + split[0]);
            return Collections.emptySet();
        }
        return Collections.singleton(new Tuple2<>(split[0], split[1]));
    }

}
