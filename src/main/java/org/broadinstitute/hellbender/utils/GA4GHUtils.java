package org.broadinstitute.hellbender.utils;

import java.nio.charset.StandardCharsets;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Base64;
import java.util.Arrays;

public class GA4GHUtils {

    public static String encodeAsSha512t24u(String input) {
        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-512");
            byte[] hashBytes = digest.digest(input.getBytes(StandardCharsets.UTF_8));
            byte[] truncatedHash = Arrays.copyOf(hashBytes, 24);
            return encodeBase64URL(truncatedHash);

        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-512 algorithm not available.", e);
        }
    }

    private static String encodeBase64URL(byte[] input) {
        return Base64.getUrlEncoder().withoutPadding().encodeToString(input);
    }

}
