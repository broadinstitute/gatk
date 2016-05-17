package org.broadinstitute.hellbender.utils.pairhmm;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.genotyper.LikelihoodMatrix;
import org.broadinstitute.hellbender.utils.haplotype.Haplotype;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Class for performing the pair HMM for local alignment using AVX instructions contained in a native shared library.
 */
public final class VectorLoglessPairHMM extends LoglessPairHMM {

    private static final Logger logger = LogManager.getLogger(VectorLoglessPairHMM.class);
    final static Boolean runningOnMac = System.getProperty("os.name", "unknown").toLowerCase().startsWith("mac");
    private static final String AVX_NATIVE_CODE_PATH_IN_JAR = "/lib/libVectorLoglessPairHMM";
    long threadLocalSetupTimeDiff = 0;
    long pairHMMSetupTime = 0;
    static Boolean isVectorLoglessPairHMMLibraryLoaded = false;

    //Hold the mapping between haplotype and index in the list of Haplotypes passed to initialize
    //Use this mapping in computeLikelihoods to find the likelihood value corresponding to a given Haplotype
    HashMap<Haplotype, Integer> haplotypeToHaplotypeListIdxMap = new HashMap<>();
    JNIHaplotypeDataHolderClass[] mHaplotypeDataArray;

    //Used to copy references to byteArrays to JNI from reads
    class JNIReadDataHolderClass {
        public byte[] readBases = null;
        public byte[] readQuals = null;
        public byte[] insertionGOP = null;
        public byte[] deletionGOP = null;
        public byte[] overallGCP = null;
    }

    //Used to copy references to byteArrays to JNI from haplotypes
    class JNIHaplotypeDataHolderClass {
        public byte[] haplotypeBases = null;
    }

    /**
     * Function to initialize the fields of JNIReadDataHolderClass and JNIHaplotypeDataHolderClass from JVM.
     * C++ codegets FieldIDs for these classes once and re-uses these IDs for the remainder of the program. Field IDs do not
     * change per JVM session
     *
     * @param readDataHolderClass      class type of JNIReadDataHolderClass
     * @param haplotypeDataHolderClass class type of JNIHaplotypeDataHolderClass
     */
    native void jniInitializeClassFields(Class<?> readDataHolderClass, Class<?> haplotypeDataHolderClass);

    /**
     * Function to report if AVX is supported on the system.
     */
    public static Boolean isAVXSupported() {
        // use a grep command to check for AVX support
        // grep exit code = 0 if a match was found
        final String command = runningOnMac ? "sysctl -a | grep machdep.cpu.features | grep -i avx" :
                                              "grep -i avx /proc/cpuinfo";
        try {
            Process process = new ProcessBuilder("/bin/sh", "-c", command).start();
            if (process.waitFor() != 0) {
                logger.warn("Error starting process to check for AVX support : " + command);
                return false;
            }
            if (process.exitValue() != 0) {
                logger.info("AVX is not supported on this system : " + command);
                return false;
            }
        }
        catch (InterruptedException | IOException e) {
            logger.warn("Error running command to check for AVX support : " + command);
            return false;
        }
        return true;
    }
    
    //The constructor is called only once inside PairHMMLikelihoodCalculationEngine
    public VectorLoglessPairHMM() throws UserException.HardwareFeatureException {
        super();

        synchronized (isVectorLoglessPairHMMLibraryLoaded) {
            // if AVX is not supported, throw an exception
            if (!isAVXSupported()) {
                throw new UserException.HardwareFeatureException("Machine does not support AVX PairHMM.");
            }

            // Load the library and initialize the FieldIDs
            if (!isVectorLoglessPairHMMLibraryLoaded) {
                try {
                    //Try loading from Java's library path first
                    //Useful if someone builds his/her own library and wants to override the bundled
                    //implementation without modifying the Java code
                    System.loadLibrary("VectorLoglessPairHMM");
                    logger.info("libVectorLoglessPairHMM found in JVM library path");
                } catch (UnsatisfiedLinkError e) {
                    //Could not load from Java's library path - try unpacking from jar
                    try {
                        logger.debug("libVectorLoglessPairHMM not found in JVM library path - trying to unpack from GATK jar file");
                        String path = AVX_NATIVE_CODE_PATH_IN_JAR + (runningOnMac ? ".dylib" : ".so");
                        File temp = IOUtils.writeTempResource(new Resource(path, VectorLoglessPairHMM.class));
                        temp.deleteOnExit();
                        System.load(temp.getAbsolutePath());
                        logger.info("libVectorLoglessPairHMM unpacked successfully from GATK jar file");
                    } catch (final Exception furtherException){
                        throw new UserException.HardwareFeatureException("Machine supports AVX, but failed to load the AVX PairHMM library from the classpath.  " +
                                "Check that your jar contains native code for your platform or choose a java HMM.", furtherException);
                    }
                }
                logger.info("Using vectorized implementation of PairHMM");
                isVectorLoglessPairHMMLibraryLoaded = true;

                //need to do this only once
                jniInitializeClassFields(JNIReadDataHolderClass.class, JNIHaplotypeDataHolderClass.class);
            }
        }
    }

    public HashMap<Haplotype, Integer> getHaplotypeToHaplotypeListIdxMap() {
        return haplotypeToHaplotypeListIdxMap;
    }

    //Used to transfer data to JNI
    //Since the haplotypes are the same for all calls to computeLikelihoods within a region, transfer the haplotypes only once to the JNI per region

    /**
     * {@inheritDoc}
     */
    @Override
    public void initialize(final List<Haplotype> haplotypes, final Map<String, List<GATKRead>> perSampleReadList,
                           final int readMaxLength, final int haplotypeMaxLength) {
        // do not need to call super.initialize()
        int numHaplotypes = haplotypes.size();
        mHaplotypeDataArray = new JNIHaplotypeDataHolderClass[numHaplotypes];
        int idx = 0;
        haplotypeToHaplotypeListIdxMap.clear();
        for (final Haplotype currHaplotype : haplotypes) {
            mHaplotypeDataArray[idx] = new JNIHaplotypeDataHolderClass();
            mHaplotypeDataArray[idx].haplotypeBases = currHaplotype.getBases();
            haplotypeToHaplotypeListIdxMap.put(currHaplotype, idx);
            ++idx;
        }
    }

    /**
     * Real compute kernel
     */
    native void jniComputeLikelihoods(int numReads, int numHaplotypes, JNIReadDataHolderClass[] readDataArray,
                                      JNIHaplotypeDataHolderClass[] haplotypeDataArray, double[] likelihoodArray, int maxNumThreadsToUse);

    /**
     * {@inheritDoc}
     */
    @Override
    public void computeLog10Likelihoods(final LikelihoodMatrix<Haplotype> logLikelihoods,
                                        final List<GATKRead> processedReads,
                                        final Map<GATKRead, byte[]> gcp) {
        if (processedReads.isEmpty()) {
            return;
        }
        if (doProfiling) {
            startTime = System.nanoTime();
        }
        int readListSize = processedReads.size();
        int numHaplotypes = logLikelihoods.numberOfAlleles();
        JNIReadDataHolderClass[] readDataArray = new JNIReadDataHolderClass[readListSize];
        int idx = 0;
        for (GATKRead read : processedReads) {
            readDataArray[idx] = new JNIReadDataHolderClass();
            readDataArray[idx].readBases = read.getBases();
            readDataArray[idx].readQuals = read.getBaseQualities();
            readDataArray[idx].insertionGOP = ReadUtils.getBaseInsertionQualities(read);
            readDataArray[idx].deletionGOP = ReadUtils.getBaseDeletionQualities(read);
            readDataArray[idx].overallGCP = gcp.get(read);
            ++idx;
        }

        mLogLikelihoodArray = new double[readListSize * numHaplotypes];      //to store results
        if (doProfiling) {
            threadLocalSetupTimeDiff = (System.nanoTime() - startTime);
        }
        //for(reads)
        //   for(haplotypes)
        //       compute_full_prob()
        jniComputeLikelihoods(readListSize, numHaplotypes, readDataArray, mHaplotypeDataArray, mLogLikelihoodArray, 12);

        int readIdx = 0;
        for (int r = 0; r < readListSize; r++) {
            int hapIdx = 0;
            for (final Haplotype haplotype : logLikelihoods.alleles()) {

                //Since the order of haplotypes in the List<Haplotype> and alleleHaplotypeMap is different,
                //get idx of current haplotype in the list and use this idx to get the right likelihoodValue
                final int idxInsideHaplotypeList = haplotypeToHaplotypeListIdxMap.get(haplotype);
                logLikelihoods.set(hapIdx, r, mLogLikelihoodArray[readIdx + idxInsideHaplotypeList]);
                ++hapIdx;
            }
            readIdx += numHaplotypes;
        }
        if (doProfiling) {
            threadLocalPairHMMComputeTimeDiff = (System.nanoTime() - startTime);
            pairHMMComputeTime += threadLocalPairHMMComputeTimeDiff;
            pairHMMSetupTime += threadLocalSetupTimeDiff;
        }
    }

    /**
     * Print final profiling information from native code
     */
    native void jniClose();

    @Override
    public void close() {
        if (doProfiling)
            logger.info("Time spent in setup for JNI call : " + (pairHMMSetupTime * 1e-9));
        super.close();
        jniClose();
    }
}
