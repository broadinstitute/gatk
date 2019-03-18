package org.broadinstitute.hellbender.utils.config;

import org.apache.spark.SparkEnv;
import org.apache.spark.TaskContext;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.TestProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.testng.annotations.Test;

import java.io.File;
import java.net.UnknownHostException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Test to make sure that the config values are honored when running on Spark.
 * Created by jonn on 4/19/18.
 */
public class SparkConfigIntegrationTest extends CommandLineProgramTest {

    //==================================================================================================================
    // Public Static Members:

    @CommandLineProgramProperties(summary = "Print the spark config.",
            oneLineSummary = "Print Config Spark Test Took",
            programGroup = TestProgramGroup.class)
    public static class PrintConfigSparkTestTool extends GATKSparkTool {

        private static final long serialVersionUID = 1L;

        static JavaSparkContext ctx;

        @Override
        protected void runTool(final JavaSparkContext ctx) {

            PrintConfigSparkTestTool.ctx = ctx;

            // Create a broadcast here:
            final Broadcast<Integer> integerBroadcast = ctx.broadcast(1);

            // Get the number of executor nodes:
            final int numExecutors = ctx.sc().statusTracker().getExecutorInfos().length;

            // Create our RDD of the size of the number of worker nodes in this run of Spark:
            final JavaRDD<Integer> javaRdd = ctx.parallelize(
                    IntStream.range(0, numExecutors).boxed().collect(Collectors.toList())
            ).repartition(numExecutors);

            // Map our RDD to our get config function so we get all the info from each worker:
            final Map<String, Properties> workerConfigMap = new HashMap<>();

            javaRdd.map(x -> getConfig( integerBroadcast ))
                    .foreach( x -> workerConfigMap.put(x.getKey(), x.getValue()) );

            // Now print each config and verify that they are all the same:
            verifyWorkerConfigsAreTheSame(workerConfigMap, numExecutors);
        }

        protected void verifyWorkerConfigsAreTheSame(final Map<String,Properties> workerConfigMap, final int numExecutors) {
            if ( workerConfigMap.size() != numExecutors ) {
                throw new GATKException("Not enough Executors were queried.  Should have been " + numExecutors + ", but was " + workerConfigMap.size() + ".");
            }

            // Get one to compare to the rest.  Which one doesn't matter because they should all be the same.
            final Map.Entry<String, Properties> firstEntry = workerConfigMap.entrySet().iterator().next();

            // Do our comparison and write out:
            final List<String> wrongWorkers = new ArrayList<>(workerConfigMap.size());
            for ( final Map.Entry<String, Properties> entry : workerConfigMap.entrySet() ) {
                System.out.println( entry.getKey() + ":" );
                entry.getValue().list(System.out);

                if ( entry.getValue() != firstEntry.getValue() ) {
                    wrongWorkers.add( entry.getKey() );
                }
            }

            if ( wrongWorkers.size() != 0 ) {
                System.out.println("The following workers had wrong entries:");
                for ( final String wrong : wrongWorkers ) {
                    System.out.print( wrong );
                }

                throw new GATKException( "Not all workers have equal properties.  These were different: " + wrongWorkers.stream().collect(Collectors.joining(",")));
            }
        }

        protected Map.Entry<String, Properties> getConfig(final Broadcast<Integer> broadcast) {

            final int partitionId = TaskContext.get().partitionId();
            final int stageId = TaskContext.get().stageId();
            final String executorId = SparkEnv.get().executorId();
            String hostAddress;
            try {
                hostAddress = java.net.InetAddress.getLocalHost().getHostAddress();
            }
            catch (final UnknownHostException ex) {
                hostAddress = "unknown";
            }

            final String workerId = hostAddress + ':' + executorId + ':' + stageId + ':' + partitionId;

            final Properties properties = new Properties();
            for ( final scala.Tuple2<String,String> confArg : ctx.getConf().getAll() ) {
                properties.put(confArg._1(), confArg._2());
            }

            return new AbstractMap.SimpleEntry<>(workerId, properties);
        }
    }

    //==================================================================================================================
    // Private Static Members:

    //==================================================================================================================
    // Private Members:

    @Test(groups = "spark")
    public void testWorkerConfigsTheSame() throws Exception {

        final File inBam = new File(getTestDataDir(), "print_reads.sorted.bam");

        final ArgumentsBuilder args = new ArgumentsBuilder();

        args.add("--" + StandardArgumentDefinitions.INPUT_LONG_NAME);
        args.add(inBam);
//        args.add("--" + StandardArgumentDefinitions.OUTPUT_LONG_NAME);
//        args.add("EMPTY_PATH");

        this.runCommandLine(args.getArgsArray());
    }

    //==================================================================================================================
    // Constructors:

    //==================================================================================================================
    // Override Methods:

    @Override
    public String getTestedClassName() {
        return PrintConfigSparkTestTool.class.getSimpleName();
    }

    //==================================================================================================================
    // Static Methods:

    //==================================================================================================================
    // Instance Methods:

    //==================================================================================================================
    // Helper Data Types:

}
