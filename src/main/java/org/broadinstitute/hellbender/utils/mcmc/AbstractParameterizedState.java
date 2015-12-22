package org.broadinstitute.hellbender.utils.mcmc;

import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Represents a collection of named {@link Parameter} objects, which may hold values of mixed type.
 * See GibbsSamplerSingleGaussianUnitTest and GibbsSamplerCopyRatioUnitTest for examples of use.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class AbstractParameterizedState {
    private final Map<String, Parameter<?>> parameterMap = new HashMap<>();

    /**
     * Constructs an AbstractParameterizedState from a List of {@link Parameter} objects with values of mixed type.
     * @param parameters    List of Parameters
     * @throws IllegalArgumentException if {@code parameters} is null or contains duplicate parameter names
     */
    public AbstractParameterizedState(final List<Parameter<?>> parameters) {
        Utils.nonNull(parameters, "List of parameters cannot be null.");
        for (final Parameter<?> parameter : parameters) {
            if (parameterMap.containsKey(parameter.name())) {
                throw new IllegalArgumentException("List of parameters cannot contain duplicate parameter names.");
            }
            parameterMap.put(parameter.name(), parameter);
        }
    }

    /**
     * Copy constructor for AbstractParameterizedState.
     * @param state AbstractParameterizedState to copy
     */
    public AbstractParameterizedState(final AbstractParameterizedState state) {
        this(new ArrayList<>(state.parameterMap.values()));
    }

    /**
     * Constructs an AbstractParameterizedState consisting of a collection of {@link Parameter} objects of identical
     * type, given a List of parameter values.  The parameter names are constructed by appending a numerical index
     * to a common parameter-name prefix.
     * @param parameterNamePrefix   common Parameter-name prefix
     * @param parameterValues       List of Parameter values
     * @param <T>                   common type for all Parameter values
     * @see     AbstractParameterizedState#initializeParameters(String, List)
     */
    public <T> AbstractParameterizedState(final String parameterNamePrefix, final List<T> parameterValues) {
        this(initializeParameters(parameterNamePrefix, parameterValues));
    }

    /**
     * Returns the value of a {@link Parameter} contained in the collection held by the AbstractParameterizedState,
     * given the parameter name and type.
     * @param parameterName         name of Parameter
     * @param parameterValueClass   class of Parameter value
     * @param <T>                   type of Parameter value
     * @return                      Parameter value
     * @throws IllegalArgumentException if {@code parameterName} does not correspond to a Parameter in the collection or
     *                                  if {@code parameterValueClass} does not match that of the Parameter in the collection
     */
    public <T> T get(final String parameterName, final Class<T> parameterValueClass) {
        try {
            return parameterValueClass.cast(parameterMap.get(parameterName).value());
        } catch (final NullPointerException | ClassCastException e) {
            if (e instanceof NullPointerException) {
                throw new IllegalArgumentException("Can only get pre-existing parameters; check parameter name.");
            }
            throw new IllegalArgumentException("Type of parameter specified in getter does not match pre-existing type.");
        }
    }

    /**
     * Abstract method for making a copy (can be shallow); must be overridden appropriately for all subclasses.
     * This can be boilerplate code; for example, for a subclass ConcreteParameterizedState, the override could read:
     * <p>{@code
     *      protected <S extends AbstractParameterizedState> S copy(final Class<S> stateClass) {
     *          return stateClass.cast(new ConcreteParameterizedState(this));
     *      }
     * }</p>
     * This copy method is primarily used by {@link ParameterizedModel} to create copies of the current state held by
     * {@link ParameterizedModel}.  These copies are ultimately stored as samples by {@link GibbsSampler}.
     * For this purpose, a shallow copy suffices if all Samplers pass new instances of parameter values to
     * {@link AbstractParameterizedState#updateParameter(String, Object)}, so that the samples
     * are not simply updated in place.  However, if this method will be used for other purposes where a deep copy
     * is required, it should be overridden accordingly.
     */
    protected abstract <S extends AbstractParameterizedState> S copy(final Class<S> stateClass);

    /**
     * Returns the set of parameter names of the Parameters contained in the collection held by the
     * AbstractParameterizedState.
     * @return  Set of Parameter names
     */
    protected Set<String> parameterNames() {
        return parameterMap.keySet();
    }

    /**
     * Updates the value of a Parameter contained in the collection held by the AbstractParameterizedState.
     * @param parameterName name of Parameter to update
     * @param value         new Parameter value
     * @param <T>           type of Parameter value
     * @throws IllegalArgumentException if {@code parameterName} does not correspond to a Parameter in the collection
     */
    protected <T> void updateParameter(final String parameterName, final T value) {
        try {
            if (!parameterMap.get(parameterName).cls().isInstance(value)) {
                throw new IllegalArgumentException("Cannot update parameter value with type different from that of current value.");
            }
            parameterMap.put(parameterName, new Parameter<>(parameterName, value));
        } catch (final NullPointerException e) {
            throw new IllegalArgumentException("Can only update pre-existing parameters; check parameter name.");
        }
    }

    /**
     * Helper method for {@link AbstractParameterizedState#AbstractParameterizedState(AbstractParameterizedState)}.
     * Returns a List of Parameters of identical type, given a List of parameter values.  The parameter names are
     * constructed by appending a numerical index to a common parameter-name prefix.
     * @param parameterNamePrefix   common Parameter-name prefix
     * @param parameterValues       List of Parameter values
     * @param <T>                   common type for all Parameter values
     * @return                      List of Parameters
     */
    private static <T> List<Parameter<?>> initializeParameters(final String parameterNamePrefix,
                                                               final List<T> parameterValues) {
        final List<Parameter<?>> initialParameters = new ArrayList<>();
        for (int i = 0; i < parameterValues.size(); i++) {
            initialParameters.add(new Parameter<>(parameterNamePrefix + i, parameterValues.get(i)));
        }
        return initialParameters;
    }
}
