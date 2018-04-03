package org.broadinstitute.hellbender.tools.spark;

import org.apache.hadoop.conf.Configuration;
import org.apache.hadoop.fs.FileSystem;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.utils.Utils;
import picard.cmdline.programgroups.OtherProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.*;
import java.nio.ByteBuffer;
import java.nio.channels.SeekableByteChannel;
import java.nio.file.Files;
import java.nio.file.NoSuchFileException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Parallel copy a file or directory from Google Cloud Storage into the HDFS file system used by Spark
 *
 * <p>This tool uses a Spark cluster to do a parallel copy of either a single file or a directory from
 * Google Cloud Storage (GCS) into the HDFS file system used by Spark to support Resilient Distributed Datasets (RDDs).
 * Files are divided into chunks of size equal to the HDFS block size (with the exception of the final
 * chunk) and each Spark task is responsible for copying one chunk. To copy all of the files in a GCS directory,
 * provide the GCS directory path, including the trailing slash. Directory copies are non-recursive so
 * subdirectories will be skipped. Within directories each file is divided into chunks independently (so this will be
 * inefficient if you have lots of files smaller than the block size). After all chunks are copied, the HDFS
 * concat method is used to stitch together chunks into single files without re-copying them.</p>
 * <p>This functionality is used by the structural variation workflow to copy reference data when a
 * Spark cluster is created, and may also be used to copy sample data to a Spark cluster.</p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li>An gcs input file or directory.</li>
 * </ul>
 *
 * <h3>Output</h3>
 * <ul>
 *     <li>A copy of the input as HDFS output.</li>
 * </ul>
 *
 * <h3>Usage example</h3>
 * <pre>
 *   gatk ParallelCopyGCSDirectoryIntoHDFSSpark \
 *     --input-gcs-path gs://bucket_name/input_reads.bam \
 *     --output-hdfs-directory hdfs://cluster-name-m:8020/directory/input_reads.bam \
 *     -- \
 *     --spark-runner GCS \
 *     --cluster my-dataproc-cluster
 * </pre>
 */
@DocumentedFeature
@CommandLineProgramProperties(
        oneLineSummary = "Parallel copy a file or directory from Google Cloud Storage into the HDFS file system used by Spark",
        summary =
        "This tool uses a Spark cluster to do a parallel copy of either a single file or a directory from" +
        " Google Cloud Storage (GCS) into the HDFS file system used by Spark to support Resilient Distributed Datasets (RDDs)." +
        " Files are divided into chunks of size equal to the HDFS block size (with the exception of the final" +
        " chunk) and each Spark task is responsible for copying one chunk. To copy all of the files in a GCS directory," +
        " provide the GCS directory path, including the trailing slash. Directory copies are non-recursive so" +
        " subdirectories will be skipped. Within directories each file is divided into chunks independently (so this will be" +
        " inefficient if you have lots of files smaller than the block size). After all chunks are copied, the HDFS" +
        " concat method is used to stitch together chunks into single files without re-copying them." +
        " This functionality is used by the structural variation workflow to copy reference data when a" +
        " Spark cluster is created, and may also be used to copy sample data to a Spark cluster.",
        programGroup = OtherProgramGroup.class)
@BetaFeature
public class ParallelCopyGCSDirectoryIntoHDFSSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    // default buffer size for reading chunks is 64MiB based on performance profiling and what appears to be conventional
    // wisdom to use a power of two for byte buffer sizes
    public static final int SIXTY_FOUR_MIB = 67108864;
    public static final String INPUT_GCS_PATH_LONG_NAME = "input-gcs-path";
    public static final String OUTPUT_HDFS_DIRECTORY_LONG_NAME = "output-hdfs-directory";
    public static final String INPUT_GLOB = "input-file-glob";
    public static final String INPUT_GLOB_ALL_FILES = "*";

    @Argument(doc = "input GCS file path (add trailing slash when specifying a directory)",
            fullName = INPUT_GCS_PATH_LONG_NAME)
    private String inputGCSPath = null;

    @Argument(doc = "optional wildcard glob to subset files in the input directory to copy",
            fullName = INPUT_GLOB)
    private String inputGlob = INPUT_GLOB_ALL_FILES;

    @Argument(doc = "output directory on HDFS to into which to transfer the data (will be created by the tool)",
            fullName = OUTPUT_HDFS_DIRECTORY_LONG_NAME)
    private String outputHDFSDirectory;

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        if (! BucketUtils.isCloudStorageUrl(inputGCSPath)) {
            throw new UserException("Input path "+ inputGCSPath + " is not a GCS URI");
        }

        if (! BucketUtils.isHadoopUrl(outputHDFSDirectory)) {
            throw new UserException("Output directory " + outputHDFSDirectory + " is not an HDFS URI");
        }

        final String inputGCSPathFinal = inputGCSPath;
        final String outputDirectoryFinal = outputHDFSDirectory;

        org.apache.hadoop.fs.Path outputHdfsDirectoryPath = new org.apache.hadoop.fs.Path(outputHDFSDirectory);

        try(FileSystem fs = outputHdfsDirectoryPath.getFileSystem(new Configuration())) {

            if (fs.exists(outputHdfsDirectoryPath)) {
                throw new UserException("Specified output directory " + outputHdfsDirectoryPath + " already exists. Please specify a new directory name.");
            }
            fs.mkdirs(outputHdfsDirectoryPath);

            final long chunkSize = getChunkSize(fs);

            final List<Path> gcsNIOPaths = getGCSFilePathsToCopy(inputGCSPathFinal, inputGlob);

            List<Tuple2<String, Integer>> chunkList = setupChunks(chunkSize, gcsNIOPaths);

            if (chunkList.size() == 0) {
                logger.info("no files found to copy");
                return;
            }

            final JavaPairRDD<String, Integer> chunkRDD = ctx.parallelizePairs(chunkList, chunkList.size());
            final JavaPairRDD<String, Tuple2<Integer, String>> chunkMappingRDD =
                    chunkRDD.mapToPair(p -> new Tuple2<>(p._1(), readChunkToHdfs(p._1(), chunkSize, p._2(), outputDirectoryFinal)));
            final Map<String, Iterable<Tuple2<Integer, String>>> chunksByFilePath = chunkMappingRDD.groupByKey().collectAsMap();

            concatenateChunks(outputDirectoryFinal, fs, gcsNIOPaths, chunksByFilePath);

        } catch (NoSuchFileException e) {
            throw new UserException("Could not locate input path " + e.getFile() + ". If you are trying to copy an entire directory, please include a trailing slash on your path.");

        } catch (IOException e) {
            throw new GATKException(e.getMessage(), e);
        }

    }

    private void concatenateChunks(final String outputDirectoryFinal, final FileSystem fs, final List<Path> gcsNIOPaths, final Map<String, Iterable<Tuple2<Integer, String>>> chunksByFilePath) throws IOException {
        for (Path path : gcsNIOPaths) {

            if (Files.isDirectory(path)) {
                continue;
            }

            final String filePath = path.toUri().toString();
            final Iterable<Tuple2<Integer, String>> chunkListForFile = chunksByFilePath.get(filePath);
            final String basename = path.getName(path.getNameCount() - 1).toString();
            final org.apache.hadoop.fs.Path outFilePath = new org.apache.hadoop.fs.Path(outputDirectoryFinal + "/" + basename);
            fs.createNewFile(outFilePath);

            SortedMap<Integer, String> chunkMap = new TreeMap<>();
            for (Tuple2<Integer, String> entry : chunkListForFile) {
                chunkMap.put(entry._1(), entry._2());
            }

            org.apache.hadoop.fs.Path[] chunkPaths = new org.apache.hadoop.fs.Path[chunkMap.size()];

            final Iterator<Integer> iterator = chunkMap.keySet().iterator();
            while (iterator.hasNext()) {
                final Integer next = iterator.next();
                final String chunkPath = chunkMap.get(next);
                chunkPaths[next] = new org.apache.hadoop.fs.Path(chunkPath);
            }

            fs.concat(outFilePath, chunkPaths);
        }
    }

    private List<Tuple2<String, Integer>> setupChunks(final long chunkSize, final List<Path> gcsNIOPaths) throws IOException {
        List<Tuple2<String, Integer>> chunkList = new ArrayList<>();
        for (Path path : gcsNIOPaths) {

            if (Files.isDirectory(path)) {
                logger.info("skipping directory " + path);
                continue;
            }

            final long fileSize = Files.size(path);
            final long chunks = fileSize / chunkSize + (fileSize % chunkSize == 0 ? 0 : 1);
            logger.info("processing path " + path + ", size = " + fileSize + ", chunks = " + chunks);

            for (int i = 0; i < chunks; i++) {
                chunkList.add(new Tuple2<>(path.toUri().toString(), i));
            }
        }
        return chunkList;
    }

    private List<Path> getGCSFilePathsToCopy(final String inputGCSPathFinal, final String inputGlob) throws IOException {
        final List<Path> gcsNIOPaths;
        final Path inputGCSNIOPath = IOUtils.getPath(inputGCSPathFinal);
        if (Files.isDirectory(inputGCSNIOPath)) {
            logger.info("transferring input directory: " + inputGCSPathFinal);
            gcsNIOPaths = Utils.stream(Files.newDirectoryStream(inputGCSNIOPath, inputGlob)).collect(Collectors.toList());
        } else {
            logger.info("transferring single file: " + inputGCSNIOPath);
            if (! INPUT_GLOB_ALL_FILES.equals(inputGlob)) {
                logger.warn("Input glob " + inputGlob + " specified, but input argument was not a directory. Ignoring glob.");
            }
            gcsNIOPaths = Collections.singletonList(inputGCSNIOPath);
        }
        return gcsNIOPaths;
    }

    static long getChunkSize(final FileSystem fs) {
        return Long.parseLong(fs.getConf().get("dfs.blocksize"));
    }

    private static final Tuple2<Integer, String> readChunkToHdfs(final String inputGCSPathFinal, final long chunkSize, final Integer chunkNum, final String outputDirectory) {
        final Path gcsPath = IOUtils.getPath(inputGCSPathFinal);
        final String basename = gcsPath.getName(gcsPath.getNameCount() - 1).toString();
        org.apache.hadoop.fs.Path outputPath = new org.apache.hadoop.fs.Path(outputDirectory);
        final String chunkPath = outputPath + "/" + basename + ".chunk." + chunkNum;

        try (SeekableByteChannel channel = Files.newByteChannel(gcsPath);
             final OutputStream outputStream = new BufferedOutputStream(BucketUtils.createFile(chunkPath))){

            final long start = chunkSize * (long) chunkNum;
            channel.position(start);
            ByteBuffer byteBuffer = ByteBuffer.allocateDirect((int) Math.min(SIXTY_FOUR_MIB, chunkSize));
            long bytesRead = 0;
            while(channel.read(byteBuffer) > 0) {
                byteBuffer.flip();
                while (byteBuffer.hasRemaining() && bytesRead < chunkSize) {
                    byte b = byteBuffer.get();
                    outputStream.write(b);
                    bytesRead++;
                }
                if (bytesRead == chunkSize) {
                    break;
                }
                if (bytesRead > chunkSize) {
                    throw new GATKException("Encountered an unknown error condition and read too many bytes; output file may be corrupt");
                }
                byteBuffer.clear();
            }
        } catch (IOException e) {
            throw new GATKException(e.getMessage() + "; inputGCSPathFinal = " + inputGCSPathFinal, e);
        }
        return new Tuple2<>(chunkNum, chunkPath);
    }

}
