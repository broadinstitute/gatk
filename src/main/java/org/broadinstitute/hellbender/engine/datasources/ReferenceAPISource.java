package org.broadinstitute.hellbender.engine.datasources;

import com.google.api.client.json.JsonFactory;
import com.google.api.services.genomics.Genomics;
import com.google.api.services.genomics.model.*;
import com.google.cloud.genomics.utils.GenomicsFactory;
import com.google.common.annotations.VisibleForTesting;
import com.google.common.primitives.Bytes;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.reference.ReferenceBases;

import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.security.GeneralSecurityException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * ReferenceAPISource makes calls to the Google Genomics API to get ReferenceBases. This is designed so it also works
 * at workers.
 * It also has helper methods for parsing Ref API URLs (intended for the Hellbender command line).
 */
public class ReferenceAPISource implements ReferenceSource, Serializable {

    // from https://cloud.google.com/genomics/data/references
    public static final String GRCH37_REF_ID = "EOSsjdnTicvzwAE";
    public static final String GRCH37_ASSEMBLY_ID ="GRCh37" ;

    // "lite" instead of the correct spelling because that's what its name is.
    public static final String GRCH37_LITE_REF_ID = "EJjur6DxjIa6KQ";
    public static final String GRCH37_LITE_ASSEMBLY_ID = "GRCh37lite";

    public static final String GRCH38_REF_ID = "EMud_c37lKPXTQ";
    public static final String GRCH38_ASSEMBLY_ID = "GRCh38";

    public static final String HG19_REF_ID = "EMWV_ZfLxrDY-wE";
    public static final String HG19_ASSEMBLY_ID = "hg19";

    public static final String HS37D5_REF_ID = "EOSt9JOVhp3jkwE";
    public static final String HS37D5_ASSEMBLY_ID = "hs37d5";

    public static final String URL_PREFIX = "gg://reference/";

    private static final long serialVersionUID = 1L;
    private static final int defaultPageSize = 1_000_000; // The number of results per request
    private static final Logger logger = LogManager.getLogger(ReferenceAPISource.class);

    // With our current design, there are three steps required (1) Create a map from reference name to Id,
    // (2) instantiate RefAPISource, and (3) call getReferenceBases, which then calls the Google Genomics API.
    // A one-off analysis showed that (1) takes ~seconds (2) ~1ms and (3) ~100ms for a small read.
    // So, creating the table is expensive and instantiating RefAPISource is cheap relative to calling getReferenceBases.
    //
    // genomicsService is created on each worker.
    private transient Genomics genomicsService;

    private Map<String, Reference> referenceMap;
    private Map<String, String> referenceNameToIdTable;

    public ReferenceAPISource( final String referenceURL ) {
        String referenceName = getReferenceSetID(referenceURL);
        this.referenceMap = getReferenceNameToReferenceTable(referenceName);
        this.referenceNameToIdTable = getReferenceNameToIdTableFromMap(referenceMap);
    }

    /**
     * Creates this ReferenceAPISource from an assembly ID by querying in the cloud APIs.
     */
    public static ReferenceAPISource fromReferenceSetAssemblyID(final String referenceSetAssemblyID){
        Utils.nonNull(referenceSetAssemblyID);
        final SearchReferenceSetsRequest content = new SearchReferenceSetsRequest();
        content.setAssemblyId(referenceSetAssemblyID);
        try {
            final Genomics genomicsService = createGenomicsService();
            final SearchReferenceSetsResponse found = genomicsService.referencesets().search(content).execute();
            final Set<String> referenceSetIds = found.getReferenceSets().stream().map(rs -> rs.getId()).collect(Collectors.toSet());
            if (referenceSetIds.isEmpty()){
                throw new UserException.UnknownReferenceSet(referenceSetAssemblyID);
            }
            if (referenceSetIds.size() > 1){
                throw new UserException.MultipleReferenceSets(referenceSetAssemblyID, referenceSetIds);
            }
            final Map<String, Reference> ret = new LinkedHashMap<>();
            for (final String rId : referenceSetIds) {
                final SearchReferencesRequest query = new SearchReferencesRequest().setReferenceSetId(rId);
                ret.putAll(genomicsService.references().search(query).execute().getReferences().stream().collect(Collectors.toMap(r -> r.getName(), r -> r)));
            }
            return new ReferenceAPISource(ret);
        } catch (final IOException e) {
            throw new UserException("Error while looking up reference set " + referenceSetAssemblyID, e);
        }
    }

    @VisibleForTesting
    Map<String, Reference> getReferenceMap() {
        return Collections.unmodifiableMap(referenceMap);
    }

    @VisibleForTesting
    public ReferenceAPISource( final Map<String, Reference> referenceMap ) {
        this.referenceMap = referenceMap;
        this.referenceNameToIdTable = getReferenceNameToIdTableFromMap(referenceMap);
    }

    public static boolean isApiSourceUrl(String url) {
        return url.startsWith(URL_PREFIX);
    }

    public static String getReferenceSetID(String url) {
        Utils.validateArg(isApiSourceUrl(url), () -> "Not a reference API source URL: "+url);
        return url.substring(URL_PREFIX.length());
    }


      /**
       * Query the Google Genomics API for reference bases spanning the specified interval from the specified
       * reference name.
       *
       * @param interval - the range of bases to retrieve.
       * @return the reference bases specified by interval and apiData (using the Google Genomics API).
       */
      @Override
      public ReferenceBases getReferenceBases(final SimpleInterval interval) {
          return getReferenceBases(interval, defaultPageSize);
      }

      /**
       * Query the Google Genomics API for reference bases spanning the specified interval from the specified
       * reference name.
       *
       * @param interval - the range of bases to retrieve.
       * @return the reference bases specified by interval and apiData (using the Google Genomics API).
       */
      public ReferenceBases getReferenceBases(final SimpleInterval interval, int pageSize) {
          Utils.nonNull(interval);

          if (genomicsService == null) {
              genomicsService = createGenomicsService();
          }
          if (!referenceNameToIdTable.containsKey(interval.getContig())) {
              throw new UserException("Contig " + interval.getContig() + " not in our set of reference names for this reference source");
          }

          try {
              final Genomics.References.Bases.List listRequest = genomicsService.references().bases().list(referenceNameToIdTable.get(interval.getContig())).setPageSize(pageSize);
              listRequest.setStart(interval.getGA4GHStart());
              listRequest.setEnd(interval.getGA4GHEnd());
              ListBasesResponse result = listRequest.execute();
              if (result.getSequence() == null) {
                  throw new UserException("No reference bases returned in query for interval " + interval + ". Is the interval valid for this reference?");
              }
              byte[] received = result.getSequence().getBytes();
              byte[] bases = received;
              if (received.length < interval.size()) {
                  final List<byte[]> blobs = new ArrayList<>();
                  blobs.add(received);
                  while (result.getNextPageToken() != null  && !result.getNextPageToken().isEmpty()) {
                      listRequest.setPageToken(result.getNextPageToken());
                      result = listRequest.execute();
                      blobs.add(result.getSequence().getBytes());
                  }
                  final byte[][] resultsArray = blobs.toArray(new byte[blobs.size()][]);
                  bases = Bytes.concat(resultsArray);
              }
              if (bases.length != interval.size()){
                  throw new UserException.ReferenceAPIReturnedUnexpectedNumberOfBytes(interval, bases);
              }
              return new ReferenceBases(bases, interval);
          } catch (IOException e) {
              throw new UserException("Query to genomics service failed for reference interval " + interval, e);
          }
      }


    /**
     * Return a sequence dictionary for the reference.
     * @param optReadSequenceDictionaryToMatch - (optional) the sequence dictionary of the reads, we'll match its order if possible.
     * @return sequence dictionary for the reference
     */
    @Override
    public SAMSequenceDictionary getReferenceSequenceDictionary(final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        return getReferenceSequenceDictionaryFromMap(referenceMap, optReadSequenceDictionaryToMatch);
    }

    /**
     * getReferenceNameToIdTableFromMap produces a table from reference name (sadly not standard names like GRCh37),
     * to IDs need for Google Genomics API calls.
     * @param referenceMap - return value from getReferenceNameToReferenceTable
     * @return returns a mapping from reference name to String ID.
     */
    public Map<String, String> getReferenceNameToIdTableFromMap(final Map<String, Reference> referenceMap) {
        Map<String, String> ret = new LinkedHashMap<>();
        for (Map.Entry<String,Reference> e : referenceMap.entrySet()) {
            ret.put(e.getKey(), e.getValue().getId());
        }
        return ret;
    }

    /**
     * Given a reference map, get a SAMSequenceDictionary that lists its contigs and their length.
     * Note that this does NOT use the referenceContig order from the API. Instead, we'll sort them
     * and try to match the read sequence dictionary order, if given.
     *
     * @param referenceMap - return value from getReferenceNameToReferenceTable
     * @param optReadSequenceDictionaryToMatch - (optional) the sequence dictionary of the reads, we'll match its order if possible.
     * @return a SAMSequenceDictionary that lists the referenceset's contigs and their length.
     */
    public SAMSequenceDictionary getReferenceSequenceDictionaryFromMap(final Map<String, Reference> referenceMap, final SAMSequenceDictionary optReadSequenceDictionaryToMatch) {
        SAMSequenceDictionary refDictionary = new SAMSequenceDictionary();
        ArrayList<SAMSequenceRecord> refContigs = new ArrayList<>();

        for (Map.Entry<String, Reference> e : referenceMap.entrySet()) {
            if (e.getKey()!=null && e.getValue().getLength()!=null) {
                refContigs.add(new SAMSequenceRecord(e.getKey(), e.getValue().getLength().intValue()));
            }
        }

        HashMap<String,Integer> indexBuilder = null;
        if (null!=optReadSequenceDictionaryToMatch) {
            indexBuilder = new LinkedHashMap<>();
            for (int i=0; i<optReadSequenceDictionaryToMatch.size(); i++) {
                final SAMSequenceRecord sequence = optReadSequenceDictionaryToMatch.getSequence(i);
                indexBuilder.put(sequence.getSequenceName(), i);
            }
        }
        final Map<String,Integer> optReadSequenceDictionaryToMatchIndex = indexBuilder;

        // GATK requires human contigs in karyotypic order: 1, 2, ..., 10, 11, ..., 20, 21, 22, X, Y with M either leading or trailing these contigs.
        // So we sort them.
        Collections.sort(refContigs, new Comparator<SAMSequenceRecord>() {
            @Override
            public int compare(SAMSequenceRecord o1, SAMSequenceRecord o2) {

                // if those are ordered in the readDictionary, then match that order
                if (null != optReadSequenceDictionaryToMatchIndex) {
                    if (optReadSequenceDictionaryToMatchIndex.containsKey(o1.getSequenceName()) && optReadSequenceDictionaryToMatchIndex.containsKey(o2.getSequenceName())) {
                        return optReadSequenceDictionaryToMatchIndex.get(o1.getSequenceName()).compareTo(optReadSequenceDictionaryToMatchIndex.get(o2.getSequenceName()));
                    }
                }

                // otherwise, order them in karyotypic order.
                int r1 = getRank(o1.getSequenceName());
                int r2 = getRank(o2.getSequenceName());
                if (r1 < r2) return -1;
                if (r2 < r1) return 1;
                return o1.getSequenceName().compareTo(o2.getSequenceName());
            }

            private int getRank(String name) {
                if (name.equalsIgnoreCase("x")) return 23;
                if (name.equalsIgnoreCase("y")) return 24;
                if (name.equalsIgnoreCase("m")) return 25;
                if (name.equalsIgnoreCase("mt")) return 25;
                StringBuilder b = new StringBuilder();
                for (char c : name.toCharArray()) {
                    if (Character.isDigit(c)) {
                        b.append(c);
                    }
                }
                String numsOnly = b.toString();
                if (numsOnly.isEmpty()) return 0;
                return Integer.parseInt(numsOnly);
            }
        });

        for (SAMSequenceRecord s : refContigs) {
            refDictionary.addSequence(s);
        }
        return refDictionary;
    }

    /**
     * getReferenceNameToReferenceTable produces a table from reference contig name
     * to the corresponding Reference object. The table is produced via a query to get the mapping. The query
     * currently takes ~10 seconds, so this table should be cached and no produced per query to get the reference bases.
     * This is clumsy and should be refactored with better Genomics API use (issue 643).
     * @param referenceSetID - the ID of the reference set to use
     * @return returns a mapping from reference name to String ID.
     */
    public Map<String, Reference> getReferenceNameToReferenceTable(final String referenceSetID) {
        Utils.nonNull(referenceSetID);
        fillGenomicsService();

        final Map<String, Reference> ret = new LinkedHashMap<>();
        try {
            final SearchReferencesRequest content = new SearchReferencesRequest();
            content.setReferenceSetId(referenceSetID);
            final SearchReferencesResponse found = genomicsService.references().search(content).execute();
            for (Reference r : found.getReferences()) {
                ret.put(r.getName(), r);
            }
        } catch (IOException e) {
            throw new UserException("Error while looking up reference set " + referenceSetID, e);
        }
        return ret;
    }

    // --------------------------------------------------------------------------------------------
    // non-public methods

    private void fillGenomicsService() {
        if (null==genomicsService) genomicsService = createGenomicsService();
    }

    private static Genomics createGenomicsService() {
        return GenomicsFactory.builder("GATK4").build().fromApplicationDefaultCredential();
    }

    // TODO: Move these to a CustomCoder. That will allow us to do something else (possibly better) for Spark.
    // TODO: See Issue #849.
    // implement methods for Java serialization, since Reference does not implement Serializable
    private void writeObject(ObjectOutputStream stream) throws IOException {
        JsonFactory jsonFactory = com.google.api.client.googleapis.util.Utils.getDefaultJsonFactory();
        stream.writeInt(referenceMap.size());
        for (Map.Entry<String, Reference> e : referenceMap.entrySet()) {
            stream.writeUTF(e.getKey());
            stream.writeUTF(jsonFactory.toString(e.getValue()));
        }
        stream.writeObject(referenceNameToIdTable);
    }
    
    private void readObject(ObjectInputStream stream) throws IOException, ClassNotFoundException {
        JsonFactory jsonFactory = com.google.api.client.googleapis.util.Utils.getDefaultJsonFactory();
        final Map<String, Reference> refs = new LinkedHashMap<>();
        int size = stream.readInt();
        for (int i = 0; i < size; i++) {
            refs.put(stream.readUTF(), jsonFactory.fromString(stream.readUTF(), Reference.class));
        }
        @SuppressWarnings("unchecked")
        final Map<String, String> refTable = (Map<String, String>) stream.readObject();
        this.referenceMap = refs;
        this.referenceNameToIdTable = refTable;
    }
}
