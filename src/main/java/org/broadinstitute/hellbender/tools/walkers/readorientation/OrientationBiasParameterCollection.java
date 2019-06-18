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
 * Container class for {@link OrientationBiasParameter} objects. The container will always have F1R2FilterConstants.NUM_KMERS objects.
 * That is, there are two duplicate entries for each reference context, since its reverse complement contains the same information.
 * We put two entries because it's simpler to query.
 */
public class OrientationBiasParameterCollection {
    final String sample;
    final private Map<String, OrientationBiasParameter> map = new HashMap<>(F1R2FilterConstants.NUM_KMERS);

    public OrientationBiasParameterCollection(final String sample){
        this.sample = sample;
        for (final String kmer : F1R2FilterConstants.ALL_KMERS){
            map.put(kmer, new OrientationBiasParameter(kmer, new double[F1R2FilterConstants.NUM_STATES], 0, 0));
        }
    }

    public String getSample() { return sample; }

    public Optional<OrientationBiasParameter> get(final String refContext){
        if (map.containsKey(refContext)){
            return Optional.of(map.get(refContext));
        } else {
            return Optional.empty();
        }
    }

    /**
     * Add an {@link OrientationBiasParameter} object to the collection. Same reference context may not be updated twice
     */
    public void set(final OrientationBiasParameter orientationBiasParameter, final ParameterType type){
        final String refContext = orientationBiasParameter.getReferenceContext();
        Utils.validate(map.get(refContext).getNumExamples() == 0,
                "updating an existing OrientationBiasParameter is not allowed. Ref context = " + refContext);
        Utils.validate(F1R2FilterConstants.CANONICAL_KMERS.contains(refContext), "set must be called on an orientationBiasParameter object with a canonical representation");

        map.put(refContext, orientationBiasParameter);

        final String revCompRefContext = SequenceUtil.reverseComplement(refContext);
        final OrientationBiasParameter revCompOrientationBiasParameter = orientationBiasParameter.getReverseComplement(type);
        map.put(revCompRefContext, revCompOrientationBiasParameter);
    }

    public void writeParameters(final File output, final ParameterType type){
        final List<OrientationBiasParameter> priors = new ArrayList<>(map.values());

        try (OrientationBiasParameter.ParameterTableWriter writer = new OrientationBiasParameter.ParameterTableWriter(
            IOUtils.fileToPath(output), sample, type)) {
            writer.writeAllRecords(priors);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while writing to %s.", output), e);
        }
    }

    /**
     * Contract: only read from the file created by {@code writeParameters}
     */
    public static OrientationBiasParameterCollection readParameters(final File input, final ParameterType type){
        final List<OrientationBiasParameter> parameters;
        final String sample;
        final String parameterType;

        try (OrientationBiasParameter.ParameterTableReader reader = new OrientationBiasParameter.ParameterTableReader(IOUtils.fileToPath(input))) {
            parameters = reader.toList();
            sample = reader.getMetadata().get(TableUtils.SAMPLE_METADATA_TAG);
            parameterType = reader.getMetadata().get(OrientationBiasParameter.PARAMETER_TYPE_TAG);

            if (parameters.size() != F1R2FilterConstants.NUM_KMERS){
                Utils.warnUser("Reading from a prior table that was not created by OrientationBiasParameterCollection::writeParameters");
            }

            if (! parameterType.equals(type.getLabel())) {
                throw new UserException(String.format("Parameter table mismatch: requested %s but the file is %s", type.getLabel(), parameterType));
            }
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception while reading from %s.", input), e);
        }

        final OrientationBiasParameterCollection orientationBiasParameterCollection = new OrientationBiasParameterCollection(sample);

        /**
         * We iterate through the canonical kmers instead of all reference contexts because otherwise we would
         * visit each reference context twice and get an error the second time, as updating an artifact prior that already
         * exists in the container class is prohibited. {@link OrientationBiasParameterCollection::set} automatically sets its reverse complement
         * for you.
         */
        for (final String refContext : F1R2FilterConstants.CANONICAL_KMERS){
            final Optional<OrientationBiasParameter> param = parameters.stream().filter(a -> a.getReferenceContext().equals(refContext)).findAny();
            if (!param.isPresent()){
                throw new UserException.BadInput( "OrientationBiasParameter object isn't present for reference context " + refContext + "in file " + input);
            }
            orientationBiasParameterCollection.set(param.get(), type);
        }
        return orientationBiasParameterCollection;
    }

    public int getNumUniqueContexts(){
        return (int) map.values().stream().filter(a -> a.getNumExamples() > 0).count()/2;
    }
}
