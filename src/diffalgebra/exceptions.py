class SymbolNameError(KeyError):
    pass

class RingMismatchError(TypeError):
    pass

class DefinitionError(KeyError):
    pass


class RingError(Exception):
    """Base exception for all ring-related errors."""
    pass

class RingValueError(RingError, ValueError):
    """Raised for invalid values in ring operations."""
    pass

class RingTypeError(RingError, TypeError):
    """Raised for type-related issues in ring operations."""
    pass

class RingWarning(Warning):
    """Base warning for non-critical ring issues."""
    pass


# Specific exceptions
class EmptySymbolError(RingValueError):
    """Raised when an empty string is provided as a symbol."""
    pass

class DuplicateSymbolError(RingValueError):
    """Raised when duplicate symbols are provided."""
    pass

class UnknownGeneratorError(RingValueError):
    """Raised when requesting a generator that doesn't exist."""
    pass

class IncompatibleRingsError(RingValueError):
    """Raised on attempting operations between elements of incompatible rings"""
    pass

class InvalidGeneratorError(RingTypeError):
    """Raised when an object is not a valid generator for an operation."""
    pass

class WrongRingError(RingValueError):
    """Raised when a generator belongs to a different ring."""
    pass

class NonIntegrableError(RingValueError):
    """Raised when attempting to integrate an expression which is not a total derivative."""

# Specific warning
class ReservedWordWarning(RingWarning):
    """Warning for using reserved words as symbols."""
    pass