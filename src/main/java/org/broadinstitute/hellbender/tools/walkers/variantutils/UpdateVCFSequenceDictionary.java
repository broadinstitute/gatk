package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.variant.utils.SAMSequenceDictionaryExtractor;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineException;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SequenceDictionaryUtils;
import picard.cmdline.programgroups.VariantManipulationProgramGroup;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;

/**
 * Updates the reference contigs in the header of the VCF format file, i.e. the reference dictionary, using the
 * dictionary from a variant, alignment, reference, or dictionary file.
 *
 * <p>This tool is designed to update the sequence dictionary in a variant file using a dictionary from another variant,
 * alignment, dictionary, or reference file. The dictionary must be valid, i.e. must contain a sequence record, for all
 * variants in the target file. The dictionary lines start with '##contig='.
 * </p>
 *
 * <p>
 *     By specifying both --replace and --disable-sequence-dictionary-validation, one can force replace an invalid
 *     sequence dictionary in a variant file with a valid sequence dictionary in another file.
 * </p>
 *
 * <h3>Usage example</h3>
 *
 * <h4>Use the contig dictionary from a BAM (SQ lines) to replace an existing dictionary in the header of a VCF.</h4>
 * <pre>
 * gatk UpdateVCFSequenceDictionary \
 *     -V cohort.vcf.gz \
 *     --source-dictionary sample.bam \
 *     --output cohort_replacedcontiglines.vcf.gz \
 *     --replace=true
 * </pre>
 *
 * <h4>Use a reference dictionary to add reference contig lines to a VCF without any.</h4>
 * <pre>
 * gatk UpdateVCFSequenceDictionary \
 *     -V resource.vcf.gz \
 *     --source-dictionary reference.dict \
 *     --output resource_newcontiglines.vcf.gz
 * </pre>
 *
 * <h4>Use the reference set to add contig lines to a VCF without any.</h4>
 * <pre>
 * gatk UpdateVCFSequenceDictionary \
 *     -V resource.vcf.gz \
 *     -R reference.fasta \
 *     --output resource_newcontiglines.vcf.gz
 * </pre>
 *
 * <p>
 *    The -O argument specifies the name of the updated file. The --source-dictionary argument specifies the
 *    input sequence dictionary. The --replace argument is optional, and forces the replacement of the dictionary
 *    if the input file already has a dictionary.
 * </p>
 *
 */
@CommandLineProgramProperties(
        summary = "Updates the sequence dictionary in a variant file using the dictionary from another variant, " +
                  "alignment, dictionary, or reference file. The dictionary must be valid for all variants in the " +
                  "target file.",
        oneLineSummary = "Updates the sequence dictionary in a variant file.",
        programGroup = VariantManipulationProgramGroup.class
)
@DocumentedFeature
public final class UpdateVCFSequenceDictionary extends VariantWalker {
    private static final Logger logger = LogManager.getLogger(UpdateVCFSequenceDictionary.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which updated variants should be written")
    public GATKPath outFile = null;

    public final static String DICTIONARY_ARGUMENT_NAME = "source-dictionary";
    @Argument(fullName=DICTIONARY_ARGUMENT_NAME,
            doc="A variant, alignment, dictionary, or reference file to use as a dictionary source " +
                "(optional if the sequence dictionary source is specified as a reference argument). The dictionary " +
                "presented must be valid (contain a sequence record) for each sequence that is referenced by any " +
                "variant in the input file.",
            optional=true)
    String dictionarySource; // optional since a reference can be provided as a dictionary source instead

    public final static String REPLACE_ARGUMENT_NAME = "replace";
    @Argument(fullName=REPLACE_ARGUMENT_NAME,
            doc="Force replacement of the dictionary if the input file already has a dictionary",
            optional=true)
    boolean replace = false; // optional since a reference can be provided as a dictionary source instead

    private VariantContextWriter vcfWriter;
    private SAMSequenceDictionary sourceDictionary;

    @Override
    public void onTraversalStart() {
        VCFHeader inputHeader = getHeaderForVariants();
        VCFHeader outputHeader = inputHeader == null ?
                new VCFHeader() :
                new VCFHeader(inputHeader.getMetaDataInInputOrder(), inputHeader.getGenotypeSamples()) ;
        getDefaultToolVCFHeaderLines().forEach(line -> outputHeader.addMetaDataLine(line));
        sourceDictionary = getBestAvailableSequenceDictionary();

        // If -replace is set, do not need to check the sequence dictionary for validity here -- it will still be
        // checked in our normal sequence dictionary validation. Warn and require opt-in via -replace if we're about to
        // clobber a valid sequence dictionary. Check the input file directly via the header rather than using the
        // engine, since it might dig one up from an index.
        if (!replace) {
            SAMSequenceDictionary oldDictionary =
                    inputHeader == null ? null : inputHeader.getSequenceDictionary();
            if (oldDictionary != null && !oldDictionary.getSequences().isEmpty())  {
                throw new CommandLineException.BadArgumentValue(
                        String.format(
                                "The input variant file %s already contains a sequence dictionary. " +
                                        "Use %s to force the dictionary to be replaced.",
                                getDrivingVariantsFeatureInput().getName(),
                                REPLACE_ARGUMENT_NAME
                        )
                );
            }

        }

        outputHeader.setSequenceDictionary(sourceDictionary);
        vcfWriter = createVCFWriter(outFile);
        vcfWriter.writeHeader(outputHeader);
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        // Validate each variant against the source dictionary manually
        SAMSequenceRecord samSeqRec = sourceDictionary.getSequence(vc.getContig());
        if (samSeqRec == null) {
            throw new CommandLineException.BadArgumentValue(
                String.format(
                    "The input variant file contains a variant (ID: \"%s\") with a reference to a sequence (\"%s\") " +
                    "that is not present in the provided dictionary",
                    vc.getID(),
                    vc.getContig()
                )
            );
        } else if (vc.getEnd() > samSeqRec.getSequenceLength()) {
            throw new CommandLineException.BadArgumentValue(
                String.format(
                    "The input variant file contains a variant (ID: \"%s\") with a reference to a sequence (\"%s\") " +
                    "that ends at a position (%d) that exceeds the length of that sequence (%d) in the provided dictionary",
                    vc.getID(),
                    vc.getContig(),
                    vc.getEnd(),
                    samSeqRec.getSequenceLength()
                )
            );
        }
        vcfWriter.add(vc);
    }

    /**
     * Close out the new variants file.
     */
    @Override
    public void closeTool() {
        if (vcfWriter != null) {
            vcfWriter.close();
        }
    }

    // Override getBestAvailableSequenceDictionary() to ensure that the new sequence dictionary specified by
    // the user (and not the sequence dictionary from the VCF we're trying to update) is used uniformly by all
    // callers. Otherwise, the wrong dictionary would be used when writing the index for the output vcf.
    //
    @Override
    public SAMSequenceDictionary getBestAvailableSequenceDictionary() {

        final SAMSequenceDictionary resultDictionary;

        final SAMSequenceDictionary masterDictionary = getMasterSequenceDictionary();
        if (dictionarySource == null) {
            if (masterDictionary != null) {
                // We'll accept the master dictionary if one was specified. Using the master dictionary
                // arg will result in sequence dictionary validation.
                logger.warn("Using the dictionary supplied via the {} argument", StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME);
                resultDictionary = masterDictionary;
            }
            else if (hasReference()) {
                    resultDictionary = getReferenceDictionary();
            } else {
                throw new CommandLineException.MissingArgument(
                        DICTIONARY_ARGUMENT_NAME, "A dictionary source file or reference file must be provided");
            }
        } else {
            if (masterDictionary != null) {
                throw new CommandLineException(String.format("Only one of %s or %s may be specified on the command line",
                        StandardArgumentDefinitions.SEQUENCE_DICTIONARY_NAME,
                        DICTIONARY_ARGUMENT_NAME));
            }
            resultDictionary = SAMSequenceDictionaryExtractor.extractDictionary(IOUtils.getPath(dictionarySource));
            if (resultDictionary == null || resultDictionary.getSequences().isEmpty()) {
                throw new CommandLineException.BadArgumentValue(
                    String.format(
                        "The specified dictionary source has an empty or invalid sequence dictionary: %s",
                        dictionarySource)
                );
            }
        }

        if( seqValidationArguments.performSequenceDictionaryValidation()
                && resultDictionary != null
                && dictionaryHasMissingLengths(resultDictionary)) {
            throw new UserException.SequenceDictionaryIsMissingContigLengths(dictionarySource, resultDictionary);
        }

        return resultDictionary;
    }

    private boolean dictionaryHasMissingLengths(final SAMSequenceDictionary resultDictionary) {
        return resultDictionary.getSequences().stream().anyMatch(s -> s.getSequenceLength() == SAMSequenceRecord.UNKNOWN_SEQUENCE_LENGTH);
    }
}
