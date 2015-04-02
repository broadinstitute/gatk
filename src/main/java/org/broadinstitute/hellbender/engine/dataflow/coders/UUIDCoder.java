package org.broadinstitute.hellbender.engine.dataflow.coders;

import com.google.cloud.dataflow.sdk.coders.DelegateCoder;
import com.google.cloud.dataflow.sdk.coders.StringUtf8Coder;

import java.util.UUID;

/**
 * Coder for UUIDs.
 *
 * This is basically StringDelegateCoder, except that we need to call fromString instead of a one arg constructor.
 */
public class UUIDCoder {

    public static final DelegateCoder<UUID, String> CODER =
            DelegateCoder.of(
                    StringUtf8Coder.of(),
                    new DelegateCoder.CodingFunction<UUID, String>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public String apply(UUID uuid) throws Exception {
                            return uuid.toString();
                        }
                    },
                    new DelegateCoder.CodingFunction<String, UUID>() {
                        private static final long serialVersionUID = 1l;

                        @Override
                        public UUID apply(String s) throws Exception {
                            return UUID.fromString(s);
                        }
                    }
            );
}
