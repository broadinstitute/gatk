package org.broadinstitute.hellbender.tools.walkers.featuremapping;

import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.FlowBasedKeyCodec;
import org.broadinstitute.hellbender.utils.read.FlowBasedRead;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.json.JSONObject;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;

import static org.broadinstitute.hellbender.utils.read.FlowBasedReadUtils.getReadFlowOrder;

public class SnvqrFeatureFactory {

    static final boolean PROVIDE_DEFAULT_FEATURES = false;
    private static final Logger logger = LogManager.getLogger(SnvqrFeatureFactory.class);

    // cycle skip cache
    private static String cycleSkipCacheFlowOrder;
    private static Map<String,Boolean> cycleSkipCache = new LinkedHashMap<>();

    // local conf
    private static String fillEdgeBaseWithValue = "A";

    // PPM-Seq constants
    static final double STRAND_RATIO_LOWER_THRESH = 0.27;
    static final double STRAND_RATIO_UPPER_THRESH = 0.73;
    static final int MIN_TOTAL_HMER_LENGTHS_IN_TAGS = 4;
    static final int MAX_TOTAL_HMER_LENGTHS_IN_TAGS = 8;
    static final int MIN_STEM_END_MATCHED_LENGTH = 11;

    private static Map<String, SnvqrFeature> features = new LinkedHashMap<>() {
        {
            // add all FlowFeatureMapper vcf attributes
            for ( final String key : FlowFeatureMapper.gatherVcfAttributes(null, 0, 0, true).keySet() ) {
                put(key, new SnvqrFeature(key) {
                    @Override
                    public Object getObjectValue(MappedFeature mappedFeature, byte nonCalledBase) {
                        if ( mappedFeature.vcfAttrCache == null )
                            mappedFeature.vcfAttrCache = FlowFeatureMapper.gatherVcfAttributes(mappedFeature, 0, 0, true);

                        // X_SCORE is a special case
                        if ( key.equals("X_SCORE") ) {
                            final String altKey = key + "_" + (char)nonCalledBase;
                            //return mappedFeature.vcfAttrCache.get(mappedFeature.vcfAttrCache.containsKey(altKey) ? altKey : key);
                            return mappedFeature.vcfAttrCache.containsKey(altKey) ? mappedFeature.vcfAttrCache.get(altKey) : 0;
                        } else {
                            return mappedFeature.vcfAttrCache.get(key);
                        }
                    }

                    @Override
                    public boolean isNonCalledBasedIndependent() {
                        return !getName().startsWith("X_SCORE");
                    }
                });
            }

            // read attributes
            for ( final String key : "rq".split(",") ) {
                put(key, new SnvqrFeature(key) {
                    @Override
                    public Object getObjectValue(MappedFeature mappedFeature) {
                        return mappedFeature.read.getAttributeAsString(key);
                    }
                    @Override
                    public boolean isMappedFeatureIndependent() {
                        return true;
                    }
                });
            }

            // max_softclip_length
            put("max_softclip_length", new SnvqrFeature("max_softclip_length") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    int front = (mappedFeature.read.getCigar().getFirstCigarElement().getOperator() == CigarOperator.SOFT_CLIP)
                                        ? mappedFeature.read.getCigar().getFirstCigarElement().getLength() : 0;
                    int back = (mappedFeature.read.getCigar().getLastCigarElement().getOperator() == CigarOperator.SOFT_CLIP)
                            ? mappedFeature.read.getCigar().getLastCigarElement().getLength() : 0;
                    return Math.max(front, back);
                }
                @Override
                public boolean isMappedFeatureIndependent() {
                    return true;
                }
            });

            // ref (which is actually the called_base)
            put("ref", new SnvqrFeature("ref") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    return String.valueOf((char)mappedFeature.readBases[0]);
                }
                @Override
                public boolean isNonCalledBasedIndependent() {
                    return true;
                }
            });

            // alt (which is actually the non_called_base)
            put("alt", new SnvqrFeature("alt") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature, byte nonCalledBase) {
                    return String.valueOf((char)nonCalledBase);
                }
            });

            // is_forward
            put("is_forward", new SnvqrFeature("is_forward") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    return true;
                }
                @Override
                public boolean isMappedFeatureIndependent() {
                    return true;
                }
            });

            // is_cycle_skip
            put("is_cycle_skip", new SnvqrFeature("is_cycle_skip") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature, byte nonCalledBase, SAMFileHeader header) {

                    // establish flow order
                    String flowOrder = new String(getReadFlowOrder(header, mappedFeature.read));

                    // calculate left+right motif
                    final byte[] bases = mappedFeature.read.getBasesNoCopy();
                    final int motifSize = 1;
                    final int leftMotifStart = Math.max(0, mappedFeature.readBasesOffset - motifSize);
                    final int rightMotifEnd = Math.min(bases.length, mappedFeature.readBasesOffset + motifSize + 1);
                    final byte[] mofitBases = Arrays.copyOfRange(bases, leftMotifStart, rightMotifEnd);

                    return isCycleSkip(mofitBases, nonCalledBase, mappedFeature.readBasesOffset - leftMotifStart, flowOrder);
                }
            });

            // next1-5, prev1-5
            int[] range = {1,2,3,4,5};
            for ( final int dist : range ) {
                put("prev" + dist, new SnvqrFeature("prev" + dist) {
                    @Override
                    public Object getObjectValue(MappedFeature mappedFeature) {
                        int ofs = mappedFeature.readBasesOffset - dist;
                        return (ofs >= 0) ? String.valueOf((char)mappedFeature.read.getBasesNoCopy()[ofs]) : fillEdgeBaseWithValue;
                    }
                    @Override
                    public boolean isNonCalledBasedIndependent() {
                        return true;
                    }
                });
                put("next" + dist, new SnvqrFeature("next" + dist) {
                    @Override
                    public Object getObjectValue(MappedFeature mappedFeature) {
                        int ofs = mappedFeature.readBasesOffset + dist;
                        return (ofs < mappedFeature.read.getBasesNoCopy().length) ? String.valueOf((char)mappedFeature.read.getBasesNoCopy()[ofs]) : fillEdgeBaseWithValue;
                    }
                    @Override
                    public boolean isNonCalledBasedIndependent() {
                        return true;
                    }
                });
            }

            // strand_ratio_category_start
            put("strand_ratio_category_start", new SnvqrFeature("strand_ratio_category_start") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    return calcStrandRatioCategory(mappedFeature.read, "as", "ts");
                }
                @Override
                public boolean isMappedFeatureIndependent() {
                    return true;
                }
            });

            // strand_ratio_category_end
            put("strand_ratio_category_end", new SnvqrFeature("strand_ratio_category_end") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    if ( !mappedFeature.read.hasAttribute("a3") ) {
                        return "END_UNREACHED";
                    } else {
                        return calcStrandRatioCategory(mappedFeature.read, "ae", "te");
                    }
                }
                @Override
                public boolean isMappedFeatureIndependent() {
                    return true;
                }
            });

            // hmer_context_ref
            put("hmer_context_ref", new SnvqrFeature("hmer_context_ref") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    return hmerContext(mappedFeature, mappedFeature.readBases[0]);
                }
                @Override
                public boolean isNonCalledBasedIndependent() {
                    return true;
                }
            });

            // hmer_context_alt
            put("hmer_context_alt", new SnvqrFeature("hmer_context_alt") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature, byte nonCalledBase) {
                    return hmerContext(mappedFeature, nonCalledBase);
                }
            });

            // X_PHRED_SCORE
            put("X_PHRED_SCORE", new SnvqrFeature("X_PHRED_SCORE") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    return 10.0 * mappedFeature.score;
                }
                @Override
                public boolean isNonCalledBasedIndependent() {
                    return true;
                }
            });

            // distance_from_edge
            put("distance_from_edge", new SnvqrFeature("distance_from_edge") {
                @Override
                public Object getObjectValue(MappedFeature mappedFeature) {
                    final int index = mappedFeature.index;
                    final int length = mappedFeature.read.getLength();
                    return Math.min(index, length - index);
                }
                @Override
                public boolean isNonCalledBasedIndependent() {
                    return true;
                }
            });
        }
    };

    private synchronized static Object isCycleSkip(byte[] mofitBases, byte nonCalledBase, int nonCalledBaseOffset, String flowOrder) {

        // construct cache key
        StringBuilder key = new StringBuilder(new String(mofitBases));
        key.append((char)nonCalledBase);
        if ( !flowOrder.equals(cycleSkipCacheFlowOrder) ) {
            cycleSkipCacheFlowOrder = flowOrder;
            cycleSkipCache.clear();
        }
        Boolean b = cycleSkipCache.get(key.toString());
        if ( b != null ) {
            return b;
        }

        // create reference and alt keys
        final int[] refKey = FlowBasedKeyCodec.baseArrayToKey(mofitBases, flowOrder);
        mofitBases[nonCalledBaseOffset] = nonCalledBase;
        final int[] altKey = FlowBasedKeyCodec.baseArrayToKey(mofitBases, flowOrder);

        // assign initial css
        boolean cssValue = (refKey.length != altKey.length);

        // if same length (NS) then see if it is possible-cycle-skip
        if ( !cssValue ) {
            for ( int n = 0 ; n < refKey.length ; n++ ) {
                if ( (refKey[n] == 0) ^ (altKey[n] == 0) ) {
                    cssValue = true;
                    break;
                }
            }
        }

        // enter into cache
        cycleSkipCache.put(key.toString(), cssValue);
        return cssValue;
    }

    private static Map<String, String> aliases = new LinkedHashMap<>() {
        {
            put("called_base", "alt");
            put("non_called_base", "ref");
            put("read_base", "alt");
            put("candidate_base", "ref");

            put("hmer_context_called_base", "hmer_context_alt");
            put("hmer_context_non_called_base", "hmer_context_ref");
            put("hmer_context_read_base", "hmer_context_alt");
            put("hmer_context_candidate_base", "hmer_context_ref");

            int[] range = {1,2,3,4,5};
            for ( final int dist : range ) {
                put("called_base_" + dist + "bp_before", "prev" + dist);
                put("called_base_" + dist + "bp_after", "next" + dist);
                put("read_base_" + dist + "bp_before", "prev" + dist);
                put("read_base_" + dist + "bp_after", "next" + dist);
            }
        }
    };

    public static SnvqrFeature getFeature(final String featureName, final JSONObject conf) {

        // translate alias?
        final String name = aliases.containsKey(featureName) ? aliases.get(featureName) : featureName;

        // find feature
        if ( !features.containsKey(name) ) {
            if ( !PROVIDE_DEFAULT_FEATURES ) {
                throw new GATKException("no such features: " + name);
            } else {
                logger.warn("providing default feature implementation for: " + name);
                return new SnvqrFeature(name) {
                    @Override
                    public Object getObjectValue(MappedFeature mappedFeature) {
                        return 0;
                    }
                };
            }
        }

        // locate feature
        final SnvqrFeature snvqrFeature = features.get(name);

        // configure category?
        if ( conf != null ) {
            if ( conf.has("categorical_features") ) {
                JSONObject categoricalFeatures = conf.getJSONObject("categorical_features");
                if (categoricalFeatures.has(name)) {
                    JSONObject dict = categoricalFeatures.getJSONObject(name);
                    snvqrFeature.setCategoryDictionary(dict.toMap());
                } else if (!name.equals(featureName) && categoricalFeatures.has(featureName)) {
                    JSONObject dict = categoricalFeatures.getJSONObject(featureName);
                    snvqrFeature.setCategoryDictionary(dict.toMap());
                }
            }
        }

        // make sure the feature is named as requested
        snvqrFeature.setName(featureName);

        return snvqrFeature;
    }

    private static String calcStrandRatioCategory(GATKRead read, String aAttrName, String tAttrName) {
        final float a = read.hasAttribute(aAttrName) ? read.getAttributeAsFloat(aAttrName) : 0;
        final float t = read.hasAttribute(tAttrName) ? read.getAttributeAsFloat(tAttrName) : 0;
        if ( (a + t) < MIN_TOTAL_HMER_LENGTHS_IN_TAGS || (a + t) > MAX_TOTAL_HMER_LENGTHS_IN_TAGS ) {
            return "UNDETERMINED";
        } else {
            if ( a == 0 ) {
                return "LIG";
            } else if ( t == 0 ) {
                return "HYB";
            } else {
                final float ratio = t / (a + t);
                if ( STRAND_RATIO_LOWER_THRESH <= ratio && ratio <= STRAND_RATIO_UPPER_THRESH ) {
                    return "MIXED";
                } else {
                    return "UNDETERMINED";
                }
            }
        }
    }

    private static int hmerContext(MappedFeature mappedFeature, byte base) {
        final byte[] bases = mappedFeature.read.getBasesNoCopy();
        int hmer = 1;
        for ( int offset = mappedFeature.readBasesOffset + 1 ; offset < bases.length && bases[offset] == base ; offset++ ) {
            hmer++;
        }
        return Math.min(hmer, FlowBasedRead.MAX_CLASS);
    }

    public static void setFillEdgeBaseWithValue(String fillEdgeBaseWithValue) {
        SnvqrFeatureFactory.fillEdgeBaseWithValue = fillEdgeBaseWithValue;
    }

}
