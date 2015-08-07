package org.broadinstitute.hellbender.engine.dataflow.datasources;

import com.google.cloud.dataflow.sdk.transforms.SerializableFunction;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.Serializable;
import java.util.Collections;
import java.util.Map;

/**
 * RefAPIMetadata contains all of the information necessary to make a call to the Google Genomics API, currently
 * it's only the reference name and the table from reference name to id.
 */
public class RefAPIMetadata implements Serializable {
    private static final long serialVersionUID = 1L;

    /**
     * @param referenceName name of the reference
     * @param referenceNameToIdTable table from reference name to id
     */
    public RefAPIMetadata(String referenceName, Map<String, String> referenceNameToIdTable) {
        this(referenceName, referenceNameToIdTable, RefWindowFunctions.IDENTITY_FUNCTION);
    }

    /**
     * @param referenceName name of the reference
     * @param referenceNameToIdTable table from reference name to id
     * @param referenceWindowFunction custom reference window function used to map reads to desired reference bases
     */
    public RefAPIMetadata(String referenceName, Map<String, String> referenceNameToIdTable, SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction) {
        this.referenceName = referenceName;
        this.referenceNameToIdTable = Collections.unmodifiableMap(referenceNameToIdTable);
        this.referenceWindowFunction = referenceWindowFunction;
    }

    private final String referenceName;
    private final Map<String, String> referenceNameToIdTable;
    private final SerializableFunction<GATKRead, SimpleInterval> referenceWindowFunction;

    public String getReferenceName() {
        return referenceName;
    }

    public Map<String, String> getReferenceNameToIdTable() {
        return referenceNameToIdTable;
    }

    public SerializableFunction<GATKRead, SimpleInterval> getReferenceWindowFunction() {
        return referenceWindowFunction;
    }
}
