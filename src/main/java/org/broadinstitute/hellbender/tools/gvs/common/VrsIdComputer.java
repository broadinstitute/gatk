package org.broadinstitute.hellbender.tools.gvs.common;

import org.broadinstitute.hellbender.exceptions.GATKException;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Base64;

/**
 * Compute VRS (Variation Representation Specification) identifiers.
 * 
 * This class implements the three-step VRS ID algorithm:
 * 1. ga4gh_serialize() - canonical JSON from inherent fields only
 * 2. sha512t24u() - SHA-512 truncated to 24 bytes, base64url-encoded
 * 3. ga4gh_identify() - prepend "ga4gh:TYPE." prefix
 * 
 * Reference: VRS 2.0.1 specification
 */
public class VrsIdComputer {

    /**
     * Compute the VRS Allele identifier for a normalized allele.
     *
     * @param refgetAccession GA4GH RefGet accession (e.g., "SQ.IIB53T8CNeJJdUqzn9V_JnRtQadwWCbl")
     * @param start          0-based interbase start coordinate
     * @param end            0-based interbase end coordinate
     * @param altSequence    Alternate allele sequence (e.g., "T" for SNP, "AAGT" for insertion)
     * @return               VRS Allele ID (e.g., "ga4gh:VA.0AePZIWZUNsUlQTamyLrjm2HWUw2opLt")
     */
    public static String computeAlleleId(
            String refgetAccession, long start, long end, String altSequence) {
        String locationDigest = computeLocationDigest(refgetAccession, start, end);
        String alleleDigest = computeAlleleDigest(locationDigest, altSequence);
        return "ga4gh:VA." + alleleDigest;
    }

    /**
     * Compute the VRS SequenceLocation identifier.
     *
     * @param refgetAccession GA4GH RefGet accession
     * @param start          0-based interbase start coordinate
     * @param end            0-based interbase end coordinate
     * @return               VRS SequenceLocation ID (e.g., "ga4gh:SL.wIlaGykfwHIpPY2Fcxtbx4TINbbODFVz")
     */
    public static String computeLocationId(
            String refgetAccession, long start, long end) {
        return "ga4gh:SL." + computeLocationDigest(refgetAccession, start, end);
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Internal digest computation
    // ──────────────────────────────────────────────────────────────────────────

    private static String computeLocationDigest(String refgetAccession, long start, long end) {
        // SequenceReference inherent fields: refgetAccession, type (sorted)
        String seqRefJson = "{\"refgetAccession\":" + jsonString(refgetAccession)
                          + ",\"type\":\"SequenceReference\"}";

        // SequenceLocation inherent fields: end, sequenceReference, start, type (sorted)
        String locationJson = "{"
                + "\"end\":" + end + ","
                + "\"sequenceReference\":" + seqRefJson + ","
                + "\"start\":" + start + ","
                + "\"type\":\"SequenceLocation\""
                + "}";

        return sha512t24u(locationJson.getBytes(StandardCharsets.UTF_8));
    }

    private static String computeAlleleDigest(String locationDigest, String altSequence) {
        // LiteralSequenceExpression inherent fields: sequence, type (sorted)
        String stateJson = "{\"sequence\":" + jsonString(altSequence)
                         + ",\"type\":\"LiteralSequenceExpression\"}";

        // Allele inherent fields: location, state, type (sorted)
        // location is Ga4ghIdentifiableObject → replaced by bare digest string
        String alleleJson = "{"
                + "\"location\":" + jsonString(locationDigest) + ","
                + "\"state\":" + stateJson + ","
                + "\"type\":\"Allele\""
                + "}";

        return sha512t24u(alleleJson.getBytes(StandardCharsets.UTF_8));
    }

    // ──────────────────────────────────────────────────────────────────────────
    // Cryptographic digest (SHA-512 truncated to 24 bytes)
    // ──────────────────────────────────────────────────────────────────────────

    /**
     * Compute sha512t24u digest: SHA-512 → leftmost 24 bytes → base64url (no padding).
     * Produces exactly 32 ASCII characters.
     *
     * @param data UTF-8 bytes to hash
     * @return 32-character base64url digest
     */
    private static String sha512t24u(byte[] data) {
        try {
            MessageDigest md = MessageDigest.getInstance("SHA-512");
            byte[] fullDigest = md.digest(data);
            byte[] truncated = new byte[24];
            System.arraycopy(fullDigest, 0, truncated, 0, 24);
            return Base64.getUrlEncoder().withoutPadding().encodeToString(truncated);
        } catch (NoSuchAlgorithmException e) {
            throw new GATKException("SHA-512 algorithm not available", e);
        }
    }

    // ──────────────────────────────────────────────────────────────────────────
    // JSON helpers
    // ──────────────────────────────────────────────────────────────────────────

    /**
     * Produce a JSON string literal with proper escaping per RFC 8785 §3.2.2.2.
     * For VRS, field values are ASCII-only (RefGet accessions, type literals, DNA sequences).
     *
     * @param value string to encode
     * @return JSON string literal (quoted and escaped)
     */
    private static String jsonString(String value) {
        StringBuilder sb = new StringBuilder("\"");
        for (char c : value.toCharArray()) {
            switch (c) {
                case '"':  sb.append("\\\""); break;
                case '\\': sb.append("\\\\"); break;
                case '\b': sb.append("\\b");  break;
                case '\f': sb.append("\\f");  break;
                case '\n': sb.append("\\n");  break;
                case '\r': sb.append("\\r");  break;
                case '\t': sb.append("\\t");  break;
                default:
                    if (c < 0x20) {
                        sb.append(String.format("\\u%04x", (int) c));
                    } else {
                        sb.append(c);
                    }
            }
        }
        sb.append("\"");
        return sb.toString();
    }
}
