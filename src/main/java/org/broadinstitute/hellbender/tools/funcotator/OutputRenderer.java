package org.broadinstitute.hellbender.tools.funcotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFHeaderLineType;
import htsjdk.variant.vcf.VCFInfoHeaderLine;
import org.broadinstitute.hellbender.tools.funcotator.dataSources.TableFuncotation;
import org.broadinstitute.hellbender.tools.funcotator.metadata.VcfFuncotationMetadata;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.LinkedHashMap;
import java.util.List;
import java.util.stream.Collectors;

/**
 * An abstract class to allow for writing output for the Funcotator.
 * Used to output Funcotations to a location of the user's choice.
 * For example, writing out to a VCF file (e.g. {@link org.broadinstitute.hellbender.tools.funcotator.vcfOutput.VcfOutputRenderer}).
 * It is also possible to write to a stream, or other connection (e.g. a UDP port)
 * at the user's discretion.
 * Created by jonn on 8/30/17.
 */
public abstract class OutputRenderer implements AutoCloseable {

    //==================================================================================================================
    /**
     * {@link LinkedHashMap} of manually specified annotations to add to each output in addition to annotations provided
     * to {@link OutputRenderer#write(VariantContext, FuncotationMap)}.
     */
    protected LinkedHashMap<String, String> manualAnnotations;

    /**
     * {@link List} of the {@link DataSourceFuncotationFactory} objects that are being used in this run of {@link Funcotator}.
     */
    protected List<DataSourceFuncotationFactory> dataSourceFactories;

    /**
     * The version of the tool used to produce the output file.
     */
    protected final String toolVersion;

    //==================================================================================================================

    /**
     * Initialize an OutputRenderer
     *
     * @param toolVersion The version number of the tool used to produce the output file (must not be null).
     */
    public OutputRenderer(final String toolVersion) {
        this.toolVersion = Utils.nonNull(toolVersion);
    }

    /**
     * @return the version number of the tool used to produce the output file
     */
    public String getToolVersion() {
        return toolVersion;
    }

    /**
     * @return A {@link String} containing information about the data sources that are used to create the {@link Funcotation}s by this {@link OutputRenderer}.
     */
    public String getDataSourceInfoString() {

        return dataSourceFactories.stream()
                .map(DataSourceFuncotationFactory::getInfoString)
                .collect(Collectors.joining(" | "));
    }

    /**
     * Close the {@link OutputRenderer}.
     */
    public abstract void close();

    /**
     * Write the given {@code variant} and {@code txToFuncotationMap} to the output file.
     * @param variant {@link VariantContext} to write to the file.
     * @param txToFuncotationMap {@link FuncotationMap} to add to the given {@code variant} on output.
     */
    public abstract void write(final VariantContext variant, final FuncotationMap txToFuncotationMap);

    /**
     * Utility for output renderers.
     *
     * Given a {@link LinkedHashMap} of field:value pairs (String:String), create a funcotation (funcotation metadata included!)
     *  with the best representation possible.
     *
     * @param data Never {@code null}
     * @param altAllele Never {@code null}
     * @param datasourceName Name to use as the datasource.  Never {@code null}
     * @return A funcotation with all fields considered strings and a generic description.  Funcotation metadata is populated.
     *   Never {@code null}
     */
    public static Funcotation createFuncotationFromLinkedHashMap(final LinkedHashMap<String, String> data, final Allele altAllele, final String datasourceName) {
        Utils.nonNull(data);
        Utils.nonNull(altAllele);
        Utils.nonNull(datasourceName);
        final List<VCFInfoHeaderLine> manualAnnotationHeaderLines = data.entrySet().stream()
                .map(e -> new VCFInfoHeaderLine(e.getKey(), 1, VCFHeaderLineType.String, "Specified from map: " + e.getKey() + ":" + e.getValue() ))
                .collect(Collectors.toList());
        return TableFuncotation.create(data, altAllele, datasourceName, VcfFuncotationMetadata.create(manualAnnotationHeaderLines));
    }
}
