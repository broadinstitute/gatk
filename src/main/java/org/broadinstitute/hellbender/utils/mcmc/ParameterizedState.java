package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Represents a mapped collection of {@link Parameter} objects, i.e., named, ordered, enumerated keys associated with
 * values of mixed type via a key -> key, value map.
 * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of use.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class ParameterizedState<T extends Enum<T> & ParameterEnum> {
    private final LinkedHashMap<T, Parameter<T, ?>> parameterMap;     //must be a Map that gives a consistently ordered keySet

    /**
     * Constructs a {@link ParameterizedState} with parameters enumerated by {@link ParameterEnum}
     * from a List of {@link Parameter} objects with values of mixed type,
     * checking that all parameters are present without duplicates.
     * @param parameters    List of Parameters
     * @throws IllegalArgumentException if {@code parameters} is null, empty, contains duplicate parameter names, or missing parameter names
     */
    public ParameterizedState(final List<Parameter<T, ?>> parameters) {
        Utils.nonNull(parameters, "List of parameters cannot be null.");
        Utils.nonEmpty(parameters, "List of parameters cannot be empty.");
        //check for duplicates in list of parameters
        final int numUniqueParameters = (int) parameters.stream().map(Parameter::getName).distinct().count();
        Utils.validateArg(numUniqueParameters == parameters.size(), "List of parameters may not contain duplicates.");
        //construct parameter map
        final LinkedHashMap<T, Parameter<T, ?>> map = new LinkedHashMap<>();
        parameters.forEach(p -> map.put(p.getName(), p));
        //check that parameter-map key set matches list of parameters
        final Class<T> keyClass = parameters.get(0).getName().getDeclaringClass();
        final Set<T> keySet = EnumSet.allOf(keyClass);
        Utils.validateArg(keySet.equals(map.keySet()), "List of parameters does not contain all parameters specified by ParameterEnum.");
        parameterMap = map;
    }

    /**
     * Copy constructor.
     * @param state state to be copied
     */
    public ParameterizedState(final ParameterizedState<T> state) {
        this(state.values());
    }

    /**
     * Returns a Set of the keys of the {@link Parameter} objects that are held internally as values of the parameter map.
     * The order of the parameters in the set is the same as in the list passed to the constructor.
     * @return a List of the {@link Parameter} objects held internally as values of the parameter map
     */
    public Set<T> keySet() {
        return Collections.unmodifiableSet(parameterMap.keySet());
    }

    /**
     * Returns a List of the {@link Parameter} objects that are held internally as values of the parameter map.
     * The order of the parameters in the list is the same as in the list passed to the constructor.
     * @return a List of the {@link Parameter} objects held internally as values of the parameter map
     */
    public List<Parameter<T, ?>> values() {
        return Collections.unmodifiableList(new ArrayList<>(parameterMap.values()));
    }

    /**
     * Returns the value of a {@link Parameter} contained in the collection held by the {@link ParameterizedState},
     * given the {@link ParameterEnum} key and type of the {@link Parameter}.
     * @param parameterKey          {@link ParameterEnum} key of {@link Parameter}
     * @param parameterValueClass   class of {@link Parameter} value
     * @param <U>                   type of {@link Parameter} value
     * @return                      {@link Parameter} value
     * @throws IllegalArgumentException if {@code parameterValueClass} does not match that of the {@link Parameter} in the collection
     */
    public <U> U get(final T parameterKey, final Class<U> parameterValueClass) {
        try {
            return parameterValueClass.cast(parameterMap.get(parameterKey).getValue());
        } catch (final ClassCastException e) {
            throw new IllegalArgumentException("Type of parameter specified in getter does not match pre-existing type.");
        }
    }

    /**
     * Updates the value of a {@link Parameter} contained in the collection held by the {@link ParameterizedState}.
     * @param parameterName {@link ParameterEnum} key of {@link Parameter} to update
     * @param value         new {@link Parameter} value
     * @param <U>           type of {@link Parameter} value
     * @throws IllegalArgumentException if {@code parameterName} does not correspond to a {@link Parameter} in the collection
     */
    protected <U> void update(final T parameterName, final U value) {
        Utils.validateArg(parameterMap.get(parameterName).getValue().getClass().isInstance(value),
                "Cannot update parameter value with type different from that of current value.");
        parameterMap.put(parameterName, new Parameter<>(parameterName, value));
    }

    @SuppressWarnings("unchecked")
    protected <S extends ParameterizedState<T>> S copy() {
        return (S) new ParameterizedState<>(values());
    }
}
