package org.broadinstitute.hellbender.utils.samples;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Class for creating a temporary in memory database of samples.
 */
public class SampleDBBuilder {
    private final PedigreeValidationType validationStrictness;
    private final SampleDB sampleDB = new SampleDB();

    private final Set<Sample> samplesFromDataSources = new LinkedHashSet<>();
    private final Set<Sample> samplesFromPedigrees = new LinkedHashSet<>();

    public SampleDBBuilder(PedigreeValidationType validationStrictness) {
        this.validationStrictness = validationStrictness;
    }

    public SampleDBBuilder addSamplesFromPedigreeFiles(final List<File> pedigreeFiles) {
        for (final File pedFile : pedigreeFiles) {
            Collection<Sample> samples = addSamplesFromPedigreeArgument(pedFile);
            samplesFromPedigrees.addAll(samples);
        }

        return this;
    }

    public SampleDBBuilder addSamplesFromPedigreeStrings(final List<String> pedigreeStrings) {
        for (final String pedString : pedigreeStrings) {
            Collection<Sample> samples = addSamplesFromPedigreeArgument(pedString);
            samplesFromPedigrees.addAll(samples);
        }

        return this;
    }

    /**
     * Parse one sample file and integrate it with samples that are already there
     * Fail quickly if we find any errors in the file
     */
    private Collection<Sample> addSamplesFromPedigreeArgument(File sampleFile) {
        final PedReader reader = new PedReader();

        try {
            return reader.parse(sampleFile, getMissingFields(sampleFile), sampleDB);
        } catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(sampleFile, e);
        }
    }

    /**
     * Integrates the collection of sample names with the samples already present
     */
    public SampleDBBuilder addSamplesFromSampleNames(final Collection<String> sampleNames) {
        Utils.nonNull(sampleNames);
        for (final String sampleName : sampleNames) {
            if (sampleDB.getSample(sampleName) == null) {
                final Sample newSample = new Sample(sampleName, null, null, null, Sex.UNKNOWN);
                sampleDB.addSample(newSample);
                samplesFromDataSources.add(newSample); // keep track of data source samples
            }
        }
        return this;
    }

    private Collection<Sample> addSamplesFromPedigreeArgument(final String string) {
        final PedReader reader = new PedReader();
        return reader.parse(string, getMissingFields(string), sampleDB);
    }

    public SampleDB getFinalSampleDB() {
        validate();
        return sampleDB;
    }

    private EnumSet<PedReader.MissingPedField> getMissingFields(final Object engineArg) {
        return EnumSet.noneOf(PedReader.MissingPedField.class);
    }

    // --------------------------------------------------------------------------------
    //
    // Validation
    //
    // --------------------------------------------------------------------------------

    private void validate() {
        validatePedigreeIDUniqueness();
        if (validationStrictness != PedigreeValidationType.SILENT) {
            // check that samples in data sources are all annotated, if anything is annotated
            if (!samplesFromPedigrees.isEmpty() && ! samplesFromDataSources.isEmpty()) {
                final Set<String> sampleNamesFromPedigrees = samplesFromPedigrees.stream().map(Sample::getID).collect(Collectors.toSet());

                for (final Sample dsSample : samplesFromDataSources)
                    if (!sampleNamesFromPedigrees.contains(dsSample.getID())) {
                        throw new UserException("Sample " + dsSample.getID()
                                + " found in data sources but not in pedigree files with STRICT pedigree validation");
                    }
            }
        }
    }

    private void validatePedigreeIDUniqueness() {
        final Set<String> pedigreeIDs = samplesFromPedigrees.stream().map(Sample::getID).collect(Collectors.toSet());
        assert pedigreeIDs.size() == samplesFromPedigrees.size() :
                "The number of sample IDs extracted from the pedigree does not equal the number of samples in the pedigree. Is a sample associated with multiple families?";
    }
}
