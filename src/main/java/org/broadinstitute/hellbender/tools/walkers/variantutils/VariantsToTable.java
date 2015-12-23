package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.*;
import java.util.function.Function;

/**
 * Extract specific fields from a VCF file to a tab-delimited table
 *
 * <p>
 * This tool is designed to extract fields from the VCF to a table format that is more convenient to work with in
 * downstream analyses.</p>
 *
 * <p>The user specifies one or more
 * fields to print with the -F NAME, each of which appears as a single column in
 * the output file, with a header named NAME, and the value of this field in the VCF
 * one per line.  NAME can be any standard VCF column (CHROM, ID, QUAL) or any annotation name
 * in the INFO field (AC=10).  In addition, there are specially supported values like
 * EVENTLENGTH (length of the event), TRANSITION (for SNPs), HET (count of het genotypes),
 * HOM-REF (count of homozygous reference genotypes), HOM-VAR (count of homozygous variant
 * genotypes), NO-CALL (count of no-call genotypes), TYPE (the type of event), VAR (count of
 * non-reference genotypes), NSAMPLES (number of samples), NCALLED (number of called samples),
 * GQ (from the genotype field; works only for a file with a single sample), and MULTI-ALLELIC
 * (is the record from a multi-allelic site).  </p>
 *
 * </p>
 *
 * <h3>Input</h3>
 * <p>
 * <ul>
 *     <li>A VCF file</li>
 *     <li>A list of -F fields to write</li>
 * </ul>
 * </p>
 *
 * <h3>Output</h3>
 * <p>
 * A tab-delimited file containing the values of the requested fields in the VCF file
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 *     ./gatk-launch \
 *     VariantsToTable \
 *     -V file.vcf \
 *     -F CHROM -F POS -F ID -F QUAL -F AC \
 *     -o results.table
 * </pre>
 * <p>would produce a file that looks like:</p>
 * <pre>
 *     CHROM    POS ID      QUAL    AC
 *     1        10  .       50      1
 *     1        20  rs10    99      10
 * </pre>
 *
 * <h3>Caveat</h3>
 * <p>If a VCF record is missing a value, then the tool by default throws an error, but the special value NA can
 * be emitted instead if requested at the command line using --allowMissingData.</p>
 */
@CommandLineProgramProperties(
        summary = "This tool is designed to extract fields from the VCF to a table format " +
                "that is more convenient to work with in " +
                "downstream analyses",
        oneLineSummary = "Extract specific fields from a VCF file to a tab-delimited table",
        programGroup = VariantProgramGroup.class
)
public final class VariantsToTable extends VariantWalker {

    static final Logger logger = LogManager.getLogger(VariantsToTable.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which results should be written (defaults to stdout)")
    private String out = null;

    /**
     * -F NAME can be any standard VCF column (CHROM, ID, QUAL) or any annotation name in the INFO field (e.g., -F AC).
     * To capture GENOTYPE (FORMAT) field values, see the -GF argument. This argument accepts any number
     * of inputs.  So -F CHROM -F POS is allowed.
     */
    @Argument(fullName="fields", shortName="F", doc="The name of each field to capture for output in the table", optional=true)
    private List<String> fieldsToTake = new ArrayList<>();

    /**
     * -GF NAME can be any annotation name in the FORMAT field (e.g., GQ, PL).
     * This argument accepts any number of inputs.  So -GF GQ -GF PL is allowed.
     */
    @Argument(fullName="genotypeFields", shortName="GF", doc="The name of each genotype field to capture for output in the table", optional=true)
    private List<String> genotypeFieldsToTake = new ArrayList<>();

    /**
     * By default this tool only emits values for records where the FILTER field is either PASS or . (unfiltered).
     * Using this flag will cause VariantsToTable to emit values regardless of the FILTER field value.
     */
    @Advanced
    @Argument(fullName="showFiltered", shortName="raw", doc="If provided, field values from filtered records will be included in the output", optional=true)
    private boolean showFiltered = false;

    /**
     * By default, records with multiple ALT alleles will comprise just one line of output; note that in general this can make your resulting file
     * unreadable/malformed for certain tools like R, as the representation of multi-allelic INFO field values are often comma-separated lists
     * of values.  Using this flag will cause multi-allelic records to be split into multiple lines of output (one for each allele in the ALT field);
     * INFO field values that are not lists are copied for each of the output records while only the appropriate entry is used for lists.
     */
    @Argument(fullName="splitMultiAllelic", shortName="SMA", doc="If provided, we will split multi-allelic records into multiple lines of output", optional=true)
    private boolean splitMultiAllelic = false;

    /**
     * By default, this tool emits one line per usable VCF record (or per allele if the -SMA flag is provided).  Using the -moltenize flag
     * will cause records to be split into multiple lines of output: one for each field provided with -F or one for each combination of sample
     * and field provided with -GF. Note that, in moltenized output, all rows will have 'site' as the value of the 'Sample' column.
     */
    @Advanced
    @Argument(fullName="moltenize", shortName="moltenize", doc="If provided, we will produce molten output", optional=true)
    private boolean moltenizeOutput = false;

    /**
     * By default, this tool throws a UserException error when it encounters a record that does not contain a value for one of the requested fields.  This
     * is generally useful when you mistype -F CHRM, so that you get a friendly warning about CHROM not being
     * found before the tool runs through 40M records.  However, in some cases you genuinely want to disable this behavior, for example to allow the use
     * of fields for which not all records have a value (e.g., AC not being calculated for filtered records, if included).  When provided, this argument
     * will cause VariantsToTable to write out NA values for missing fields instead of throwing an error.
     * Note that this flag only applies to standard columns (CHROM, ID, QUAL) and the INFO field and it does not apply to the genotype field.
     */
    @Advanced
    @Argument(fullName="allowMissingData", shortName="AMD", doc="If provided, we will not require every record to contain every field", optional=true)
    private boolean allowMissingData = false;
    private static final String MISSING_DATA = "NA";


    private SortedSet<String> samples;
    private long nRecords = 0L;
    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        outputStream = createPrintStream();

        if (genotypeFieldsToTake.isEmpty()) {
            samples = Collections.emptySortedSet();
        } else {
            final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
            samples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

            // if there are no samples, we don't have to worry about any genotype fields
            if (samples.isEmpty()) {
                genotypeFieldsToTake.clear();
                logger.warn("There are no samples - the genotype fields will be ignored");
                if (fieldsToTake.isEmpty()){
                    throw new UserException("There are no samples and no fields - no output will be produced");
                }
            }
        }

        // print out the header
        if ( moltenizeOutput ) {
            outputStream.println("RecordID\tSample\tVariable\tValue");
        } else {
            final String baseHeader = Utils.join("\t", fieldsToTake);
            final String genotypeHeader = createGenotypeHeader();
            final String separator = (!baseHeader.isEmpty() && !genotypeHeader.isEmpty()) ? "\t" : "";
            outputStream.println(baseHeader + separator + genotypeHeader);
        }
    }

    private PrintStream createPrintStream() {
        try {
            return out != null ? new PrintStream(out) : System.out;
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(out, e);
        }
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if ( showFiltered || vc.isNotFiltered() ) {
            nRecords++;
            final List<List<String>> records = extractFields(vc);
            if (moltenizeOutput){
                records.forEach(record -> emitMoltenizedOutput(record));
            } else {
                records.forEach(record -> outputStream.println(Utils.join("\t", record)));
            }
        }
    }

    private static boolean isWildCard(final String s) {
        return s.endsWith("*");
    }

    private String createGenotypeHeader() {
        boolean firstEntry = true;

        final StringBuilder sb = new StringBuilder();
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                if ( firstEntry ) {
                    firstEntry = false;
                } else {
                    sb.append("\t");
                }
                // spaces in sample names are legal but wreak havoc in R data frames
                sb.append(sample.replace(" ","_"));
                sb.append('.');
                sb.append(gf);
            }
        }
        return sb.toString();
    }

    private void emitMoltenizedOutput(final List<String> record) {
        int index = 0;
        for ( final String field : fieldsToTake ) {
            outputStream.println(String.format("%d\tsite\t%s\t%s", nRecords, field, record.get(index++)));
        }
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                outputStream.println(String.format("%d\t%s\t%s\t%s", nRecords, sample.replace(" ","_"), gf, record.get(index++)));
            }
        }
    }

    /**
     * Utility function that returns the list of values for each field in fields from vc.
     *
     * @param vc                the VariantContext whose field values we can to capture
     * @return List of lists of field values
     */
    private List<List<String>> extractFields(final VariantContext vc) {

        final int numRecordsToProduce = splitMultiAllelic ? vc.getAlternateAlleles().size() : 1;
        final List<List<String>> records = new ArrayList<>(numRecordsToProduce);

        final int numFields;
        final boolean addGenotypeFields = genotypeFieldsToTake != null && !genotypeFieldsToTake.isEmpty();
        if ( addGenotypeFields ) {
            numFields = fieldsToTake.size() + genotypeFieldsToTake.size() * samples.size();
        } else {
            numFields = fieldsToTake.size();
        }

        for ( int i = 0; i < numRecordsToProduce; i++ ) {
            records.add(new ArrayList<>(numFields));
        }

        for ( final String field : fieldsToTake ) {
            if ( splitMultiAllelic && field.equals("ALT") ) { // we need to special case the ALT field when splitting out multi-allelic records
                addFieldValue(splitAltAlleles(vc), records);
            } else if ( getters.containsKey(field) ) {
                addFieldValue(getters.get(field).apply(vc), records);
            } else if ( vc.hasAttribute(field) ) {
                addFieldValue(vc.getAttribute(field, null), records);
            } else if ( isWildCard(field) ) {
                final SortedSet<String> wildVals = new TreeSet<>();
                for ( final Map.Entry<String,Object> elt : vc.getAttributes().entrySet()) {
                    if ( elt.getKey().startsWith(field.substring(0, field.length() - 1)) ) {
                        wildVals.add(elt.getValue().toString());
                    }
                }

                final String val = wildVals.isEmpty() ? MISSING_DATA : Utils.join(",", wildVals);

                addFieldValue(val, records);
            } else if ( ! allowMissingData ) {
                throw new UserException(String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc));
            } else {
                addFieldValue(MISSING_DATA, records);
            }
        }

        if ( addGenotypeFields ) {
            addGenotypeFieldsToRecords(vc, records);
        }

        return records;
    }

    private void addGenotypeFieldsToRecords(final VariantContext vc, final List<List<String>> records) {
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                if ( vc.hasGenotype(sample) && vc.getGenotype(sample).hasAnyAttribute(gf) ) {
                    if (VCFConstants.GENOTYPE_KEY.equals(gf)) {
                        addFieldValue(vc.getGenotype(sample).getGenotypeString(true), records);
                    } else {
                        addFieldValue(vc.getGenotype(sample).getAnyAttribute(gf), records);
                    }
                } else {
                    addFieldValue(MISSING_DATA, records);
                }
            }
        }
    }

    private static void addFieldValue(final Object val, final List<List<String>> result) {
        final int numResultRecords = result.size();

        // if we're trying to create a single output record, add it
        if ( numResultRecords == 1 ) {
            result.get(0).add(prettyPrintObject(val));
        }
        // if this field is a list of the proper size, add the appropriate entry to each record
        else if ( (val instanceof List) && ((List)val).size() == numResultRecords ) {
            final List<?> list = (List<?>)val;
            for ( int i = 0; i < numResultRecords; i++ ) {
                result.get(i).add(list.get(i).toString());
            }
        }
        // otherwise, add the original value to all of the records
        else {
            final String valStr = prettyPrintObject(val);
            for ( final List<String> record : result ) {
                record.add(valStr);
            }
        }
    }

    private static String prettyPrintObject(final Object val) {
        if ( val instanceof List ) {
            return prettyPrintObject(((List) val).toArray());
        }

        if ( !val.getClass().isArray() ) {
            return val.toString();
        }

        final int length = Array.getLength(val);
        if ( length == 0 ) {
            return "";
        }

        final StringBuilder sb = new StringBuilder(prettyPrintObject(Array.get(val, 0)));
        for ( int i = 1; i < length; i++ ) {
            sb.append(",");
            sb.append(prettyPrintObject(Array.get(val, i)));
        }
        return sb.toString();
    }

    // ----------------------------------------------------------------------------------------------------
    //
    //  Getting values from VC by name.
    //
    // ----------------------------------------------------------------------------------------------------

    private final Map<String, Function<VariantContext, String>> getters = new HashMap<>();
    {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", vc -> vc.getContig());
        getters.put("POS", vc -> Integer.toString(vc.getStart()));
        getters.put("REF", vc -> vc.getReference().getDisplayString());
        getters.put("ALT", vc -> {
            final StringBuilder x = new StringBuilder();
            final int n = vc.getAlternateAlleles().size();
            if ( n == 0 ) {
                return ".";
            }

            for ( int i = 0; i < n; i++ ) {
                if ( i != 0 ) {
                    x.append(",");
                }
                x.append(vc.getAlternateAllele(i));
            }
            return x.toString();
        });
        getters.put("EVENTLENGTH", vc -> {
            int maxLength = 0;
            for ( final Allele a : vc.getAlternateAlleles() ) {
                final int length = a.length() - vc.getReference().length();
                if( Math.abs(length) > Math.abs(maxLength) ) { maxLength = length; }
            }
            return Integer.toString(maxLength);
        });
        getters.put("QUAL", vc -> Double.toString(vc.getPhredScaledQual()));
        getters.put("TRANSITION", vc -> {
            if ( vc.isSNP() && vc.isBiallelic() ) {
                return GATKVariantContextUtils.isTransition(vc) ? "1" : "0";
            } else {
                return "-1";
            }
        });
        getters.put("FILTER", vc -> vc.isNotFiltered() ? "PASS" : Utils.join(",", vc.getFilters()));
        getters.put("ID", vc -> vc.getID());
        getters.put("HET", vc -> Integer.toString(vc.getHetCount()));
        getters.put("HOM-REF", vc -> Integer.toString(vc.getHomRefCount()));
        getters.put("HOM-VAR", vc -> Integer.toString(vc.getHomVarCount()));
        getters.put("NO-CALL", vc -> Integer.toString(vc.getNoCallCount()));
        getters.put("TYPE", vc -> vc.getType().toString());
        getters.put("VAR", vc -> Integer.toString(vc.getHetCount() + vc.getHomVarCount()));
        getters.put("NSAMPLES", vc -> Integer.toString(vc.getNSamples()));
        getters.put("NCALLED", vc -> Integer.toString(vc.getNSamples() - vc.getNoCallCount()));
        getters.put("MULTI-ALLELIC", vc -> Boolean.toString(vc.getAlternateAlleles().size() > 1));
    }

    private static Object splitAltAlleles(final VariantContext vc) {
        final int numAltAlleles = vc.getAlternateAlleles().size();
        if ( numAltAlleles == 1 ) {
            return vc.getAlternateAllele(0);
        }

        return vc.getAlternateAlleles();
    }
}
