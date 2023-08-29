package org.broadinstitute.hellbender.utils;

import java.io.Serializable;
import java.util.function.Supplier;

public interface SerializableSupplier <T> extends Supplier<T>, Serializable {}