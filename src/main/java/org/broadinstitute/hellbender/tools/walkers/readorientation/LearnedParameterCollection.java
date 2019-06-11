package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.*;

/**
 * Container class for {@link LearnedParameter} objects. The container will always have F1R2FilterConstants.NUM_KMERS objects.
 * That is, there are two duplicate entries for each reference context, since its reverse complement contains the same information.
 * We put two entries because it's simpler to query.
 */
public class LearnedParameterCollection {
    final String sample;
    final private Map<String, LearnedParameter> map = new HashMap<>(F1R2FilterConstants.NUM_KMERS);

    public LearnedParameterCollection(final String sample){
        this.sample = sample;
        for (final String kmer : F1R2FilterConstants.ALL_KMERS){
            map.put(kmer, new LearnedParameter(kmer, new double[F1R2FilterConstants.NUM_STATES], 0, 0));
        }
    }

    public String getSample() { return sample; }

    public Optional<LearnedParameter> get(final String refContext){
        if (map.containsKey(refContext)){
            return Optional.of(map.get(refContext));
        } else {
            return Optional.empty();
        }
    }

    /**
     * Add an artifact prior to the collection. Same reference context may not be updated twice
     */
    public void set(final LearnedParameter learnedParameter){
        final String refContext = learnedParameter.getReferenceContext();
        Utils.validate(map.get(refContext).getNumExamples() == 0,
                "updating an existing LearnedParameter is not allowed. Ref context = " + refContext);
        Utils.validate(F1R2FilterConstants.CANONICAL_KMERS.contains(refContext), "set must be called on an learnedParameter object with a canonical representation");

        map.put(refContext, learnedParameter);

        final String revCompRefContext = SequenceUtil.reverseComplement(refContext);
        final LearnedParameter revCompLearnedParameter = learnedParameter.getReverseComplement();
        map.put(revCompRefContext, revCompLearnedParameter);
    }

    public void writeArtifactPriors(final File output, ParameterType type){
        final List<LearnedParameter> priors = new ArrayList<>(map.values());

        try (LearnedParameter.LearnedParameterTableWriter writer = new LearnedParameter.LearnedParameterTableWriter(
            IOUtils.fileToPath(output), sample, type)) {
            writer.writeAllRecords(priors);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", output), e);
        }
    }

    /**
     * Contract: only read from the file created by {@code writeArtifactPriors}
     */
    public static LearnedParameterCollection readArtifactPriors(final File input, ParameterType type){
        final List<LearnedParameter> priors;
        final String sample;
        final String parameterType;

        try (LearnedParameter.LearnedParameterTableReader reader = new LearnedParameter.LearnedParameterTableReader(IOUtils.fileToPath(input))) {
            priors = reader.toList();
            sample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
            parameterType = reader.getMetadata().get(LearnedParameter.PARAMETER_TYPE_TAG);

            if (priors.size() != F1R2FilterConstants.NUM_KMERS){
                Utils.warnUser("Reading from a prior table that was not created by LearnedParameterCollection::writeArtifactPriors");
            }

            if (! parameterType.equals(type.getLabel())) {
                throw new UserException(String.format("Parameter table mismatch: requested %s but the file is %s", type.getLabel(), parameterType));
            }
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", input), e);
        }

        final LearnedParameterCollection learnedParameterCollection = new LearnedParameterCollection(sample);

        /**
         * We iterate through the canonical kmers instead of all reference contexts because otherwise we would
         * visit each reference context twice and get an error the second time, as updating an artifact prior that already
         * exists in the container class is prohibited. {@link LearnedParameterCollection.set} automatically sets its reverse complement
         * for you.
         */
        for (final String refContext : F1R2FilterConstants.CANONICAL_KMERS){
            final Optional<LearnedParameter> ap = priors.stream().filter(a -> a.getReferenceContext().equals(refContext)).findAny();
            if (!ap.isPresent()){
                throw new UserException.BadInput( "LearnedParameter object isn't present for reference context " + refContext + "in file " + input);
            }
            learnedParameterCollection.set(ap.get());
        }
        return learnedParameterCollection;
    }

    public int getNumUniqueContexts(){
        return (int) map.values().stream().filter(a -> a.getNumExamples() > 0).count()/2;
    }
}
