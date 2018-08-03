package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Container class for {@link ArtifactPrior} objects. The container will always have {@link F1R2FilterConstants.NUM_KMERS} objects.
 * That is, there are two duplicate entries for each reference context, since its reverse complement contains the same information.
 * We put two entries because it's simpler to query.
 */
public class ArtifactPriorCollection {
    final private Map<String, ArtifactPrior> map = new HashMap<>(F1R2FilterConstants.NUM_KMERS);

    public ArtifactPriorCollection(){
        for (final String kmer : F1R2FilterConstants.ALL_KMERS){
            map.put(kmer, new ArtifactPrior(kmer, new double[F1R2FilterConstants.NUM_STATES], 0, 0));
        }
    }

    public Optional<ArtifactPrior> get(final String refContext){
        if (map.containsKey(refContext)){
            return Optional.of(map.get(refContext));
        } else {
            return Optional.empty();
        }
    }

    /**
     * Add an artifact prior to the collection. Same reference context may not be updated twice
     */
    public void set(final ArtifactPrior artifactPrior){
        final String refContext = artifactPrior.getReferenceContext();
        Utils.validate(map.get(refContext).getNumExamples() == 0,
                "updating an existing ArtifactPrior is not allowed. Ref context = " + refContext);
        Utils.validate(F1R2FilterConstants.CANONICAL_KMERS.contains(refContext), "set must be called on an artifactPrior object with a canonical representation");

        map.put(refContext, artifactPrior);

        final String revCompRefContext = SequenceUtil.reverseComplement(refContext);
        final ArtifactPrior revCompArtifactPrior = artifactPrior.getReverseComplement();
        map.put(revCompRefContext, revCompArtifactPrior);
    }

    public void writeArtifactPriors(final File output){
        final List<ArtifactPrior> priors = new ArrayList<>(map.values());

        try (ArtifactPrior.ArtifactPriorTableWriter writer = new ArtifactPrior.ArtifactPriorTableWriter(output)) {
            writer.writeAllRecords(priors);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", output), e);
        }
    }

    /**
     * Contract: only read from the file created by {@code writeArtifactPriors}
     */
    public static ArtifactPriorCollection readArtifactPriors(final File input){
        final List<ArtifactPrior> priors;
        try (ArtifactPrior.ArtifactPriorTableReader reader = new ArtifactPrior.ArtifactPriorTableReader(input)) {
            priors = reader.toList();
            if (priors.size() != F1R2FilterConstants.NUM_KMERS){
                Utils.warnUser("Reading from a prior table that was not created by ArtifactPriorCollection::writeArtifactPriors");
            }
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", input), e);
        }

        final ArtifactPriorCollection artifactPriorCollection = new ArtifactPriorCollection();

        /**
         * We iterate through the canonical kmers instead of all reference contexts because otherwise we would
         * visit each reference context twice and get an error the second time, as updating an artifact prior that already
         * exists in the container class is prohibited. {@link ArtifactPriorCollection.set} automatically sets its reverse complement
         * for you.
         */
        for (final String refContext : F1R2FilterConstants.CANONICAL_KMERS){
            final Optional<ArtifactPrior> ap = priors.stream().filter(a -> a.getReferenceContext().equals(refContext)).findAny();
            if (!ap.isPresent()){
                throw new UserException.BadInput( "ArtifactPrior object isn't present for reference context " + refContext + "in file " + input);
            }
            artifactPriorCollection.set(ap.get());
        }
        return artifactPriorCollection;
    }

    public int getNumUniqueContexts(){
        return (int) map.values().stream().filter(a -> a.getNumExamples() > 0).count()/2;
    }
}
