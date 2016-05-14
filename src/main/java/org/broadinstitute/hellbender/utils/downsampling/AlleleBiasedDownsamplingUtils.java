package org.broadinstitute.hellbender.utils.downsampling;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Sets;
import htsjdk.variant.variantcontext.Allele;
import org.apache.commons.collections4.map.DefaultedMap;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.PileupElement;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.text.XReadLines;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * The purpose of this set of utilities is to downsample a set of reads to remove contamination.
 * It is NOT a general-purpose downsampling mechanism.
 */
public final class AlleleBiasedDownsamplingUtils {

    private AlleleBiasedDownsamplingUtils() {
    }

    /**
     * Computes an allele biased version of the allele counts for a given pileup
     *
     * @param alleleCounts              the allele counts for the original pileup
     * @param numReadsToRemove          number of total reads to remove per allele
     * @return non-null array of new counts needed per allele
     */
    @VisibleForTesting
    static int[] runSmartDownsampling(final int[] alleleCounts, final int numReadsToRemove) {
        final int numAlleles = alleleCounts.length;

        int maxScore = scoreAlleleCounts(alleleCounts);
        int[] alleleCountsOfMax = alleleCounts;

        final int numReadsToRemovePerAllele = numReadsToRemove / 2;

        for ( int i = 0; i < numAlleles; i++ ) {
            for ( int j = i; j < numAlleles; j++ ) {
                final int[] newCounts = alleleCounts.clone();

                // split these cases so we don't lose on the floor (since we divided by 2)
                if ( i == j ) {
                    newCounts[i] = Math.max(0, newCounts[i] - numReadsToRemove);
                } else {
                    newCounts[i] = Math.max(0, newCounts[i] - numReadsToRemovePerAllele);
                    newCounts[j] = Math.max(0, newCounts[j] - numReadsToRemovePerAllele);
                }

                final int score = scoreAlleleCounts(newCounts);

                if ( score < maxScore ) {
                    maxScore = score;
                    alleleCountsOfMax = newCounts;
                }
            }
        }

        return alleleCountsOfMax;
    }

    private static int scoreAlleleCounts(final int[] alleleCounts) {
        if ( alleleCounts.length < 2 ) {
            return 0;
        }

        final int[] alleleCountsCopy = alleleCounts.clone();
        Arrays.sort(alleleCountsCopy);

        final int maxCount = alleleCountsCopy[alleleCounts.length - 1];
        final int nextBestCount = alleleCountsCopy[alleleCounts.length - 2];
        final int remainderCount = (int)MathUtils.sum(alleleCountsCopy) - maxCount - nextBestCount;

        // try to get the best score:
        //    - in the het case the counts should be equal with nothing else
        //    - in the hom case the non-max should be zero
        return Math.min(maxCount - nextBestCount + remainderCount, Math.abs(nextBestCount + remainderCount));
    }

   /**
     *
     * Computes reads to remove based on an allele biased down-sampling
     *
     * @param alleleReadMap             original list of records per allele
     * @param contaminationFraction      the fraction of total reads to remove per allele
     * @return list of reads TO REMOVE from allele biased down-sampling
     */
    public static <A extends Allele> List<GATKRead> selectAlleleBiasedReads(final Map<A, List<GATKRead>> alleleReadMap, final double contaminationFraction) {
        Utils.nonNull(alleleReadMap, "alleleReadMap is null");
        if (contaminationFraction < 0.0 || contaminationFraction > 1.0) {
            throw new IllegalArgumentException("invalid contamination fraction " + contaminationFraction);
        }
        return selectAlleleBiasedReads(alleleReadMap, totalReads(alleleReadMap), contaminationFraction);
    }

    /**
     * Computes reads to remove based on an allele biased down-sampling
     *
     * @param alleleReadMap             original list of records per allele
     * @param contaminationFraction      the fraction of total reads to remove per allele
     * @return list of reads TO REMOVE from allele biased down-sampling
     */
    public static <A extends Allele> List<GATKRead> selectAlleleBiasedReads(final Map<A, List<GATKRead>> alleleReadMap, final int totalReads, final double contaminationFraction) {
        //no checks here - done on the public level
        final int numReadsToRemove = (int)(totalReads * contaminationFraction);

        // make a listing of allele counts
        final List<Allele> alleles = new ArrayList<>(alleleReadMap.keySet());
        alleles.remove(Allele.NO_CALL);    // ignore the no-call bin
        final int numAlleles = alleles.size();

        final int[] alleleCounts = new int[numAlleles];
        for ( int i = 0; i < numAlleles; i++ ) {
            alleleCounts[i] = alleleReadMap.get(alleles.get(i)).size();
        }

        final int[] targetAlleleCounts = runSmartDownsampling(alleleCounts, numReadsToRemove);

        final List<GATKRead> readsToRemove = new ArrayList<>(numReadsToRemove);
        for ( int i = 0; i < numAlleles; i++ ) {
            if ( alleleCounts[i] > targetAlleleCounts[i] ) {
                readsToRemove.addAll(downsampleElements(alleleReadMap.get(alleles.get(i)), alleleCounts[i] - targetAlleleCounts[i]));
            }
        }

        return readsToRemove;
    }

    /**
     * Returns the sum of length of the lists.
     */
    public static int totalReads(final Map<?, List<GATKRead>> alleleReadMap){
        return Utils.nonNull(alleleReadMap).values().stream().mapToInt(list -> list.size()).sum();
    }

    /**
     * Performs allele biased down-sampling on a pileup and computes the list of elements to remove
     *
     * @param reads                     original list of records
     * @param numElementsToRemove       the number of records to remove
     * @return the list of pileup elements TO REMOVE. The list is unmodifable.
     */
    private static List<GATKRead> downsampleElements(final List<GATKRead> reads, final int numElementsToRemove) {
        if ( numElementsToRemove == 0 ) {  //remove none
            return Collections.emptyList();
        }
        if ( numElementsToRemove >= reads.size()) {    //remove all
            return Collections.unmodifiableList(reads);
        }

        final List<GATKRead> elementsToRemove = new ArrayList<>(numElementsToRemove);
        for (final int idx : MathUtils.sampleIndicesWithoutReplacement(reads.size(), numElementsToRemove)){
            elementsToRemove.add(reads.get(idx));
        }
        return Collections.unmodifiableList(elementsToRemove);
    }

    /**
     * Create sample-contamination maps from file.
     * The format is: tab-separated with no header,
     * each line is: sampleID contaminationFraction
     *
     * @param file   Filename containing two columns: SampleID and Contamination
     * @param defaultContaminationFraction default contamination fraction, used for samples that do no specify one
     * @param sampleIDs          Set of Samples of interest (no reason to include every sample in file) or null to turn off checking
     * @param logger                      for logging output
     * @return sample-contamination Map. The returned map is a {@link DefaultedMap} that defaults to the defaultContaminationFraction for unspecified samples
     * @throws UserException if there's an IO problem reading the file.
     * @throws UserException if the file is malformed
     */

    public static DefaultedMap<String, Double> loadContaminationFile(final File file, final double defaultContaminationFraction, final Set<String> sampleIDs, final Logger logger) {
        final DefaultedMap<String, Double> sampleContamination = new DefaultedMap<>(defaultContaminationFraction);
        final Set<String> nonSamplesInContaminationFile = new LinkedHashSet<>(sampleContamination.keySet());
        try ( final XReadLines reader = new XReadLines(file, true) ){
            for (final String line : reader) {
                if (line.isEmpty()) {
                    continue;
                }

                final String [] fields = line.split("\t");
                if (fields.length != 2){
                    throw new UserException.MalformedFile("Contamination file must have exactly two, tab-delimited columns. Offending line:\n" + line);
                }
                if (fields[0].isEmpty() || fields[1].isEmpty()) {
                    throw new UserException.MalformedFile("Contamination file can not have empty strings in either column. Offending line:\n" + line);
                }

                final double contamination;
                try {
                    contamination = Double.parseDouble(fields[1]);
                } catch (final NumberFormatException e) {
                    throw new UserException.MalformedFile("Contamination file contains unparsable double in the second field. Offending line: " + line);
                }
                final String sampleName= fields[0];
                if (sampleContamination.containsKey(sampleName)) {
                    throw new UserException.MalformedFile("Contamination file contains duplicate entries for input name " + sampleName);
                }
                if (contamination < 0.0 || contamination > 1.0){
                    throw new UserException.MalformedFile("Contamination file contains unacceptable contamination value (must be 0<=x<=1): " + line);
                }
                if (sampleIDs == null || sampleIDs.contains(sampleName)) {
                    sampleContamination.put(sampleName, contamination);
                } else {
                    nonSamplesInContaminationFile.add(sampleName);
                }
            }

            //output to the user info lines telling which samples are in the Contamination File
            if (! sampleContamination.isEmpty()) {
                logger.info(String.format("The following samples were found in the Contamination file and will be processed at the contamination level therein: %s", sampleContamination.keySet().toString()));

                //output to the user info lines telling which samples are NOT in the Contamination File
                if(sampleIDs!=null){
                    final Set<String> samplesNotInContaminationFile = Sets.difference(sampleIDs, sampleContamination.keySet());
                    if (! samplesNotInContaminationFile.isEmpty()) {
                        logger.info(String.format("The following samples were NOT found in the Contamination file and will be processed at the default contamination level: %s", samplesNotInContaminationFile.toString()));
                    }
                }
            }

            //output to the user Samples that do not have lines in the Contamination File
            if (! nonSamplesInContaminationFile.isEmpty()) {
                logger.info(String.format("The following entries were found in the Contamination file but were not SAMPLEIDs. They will be ignored: %s", nonSamplesInContaminationFile.toString()));
            }

            return sampleContamination;

        } catch (IOException e) {
            throw new UserException.CouldNotReadInputFile("I/O Error while reading sample-contamination file " + file.getAbsolutePath() + ": " + e.getMessage());
        }
    }
}
