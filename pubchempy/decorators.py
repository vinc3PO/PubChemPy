import functools
import warnings
from .errors import PubChemPyDeprecationWarning


def memoized_property(fget):
    """Decorator to create memoized properties.

    Used to cache :class:`~pubchempy.Compound` and :class:`~pubchempy.Substance` properties that require an additional
    request.
    """
    attr_name = '_{0}'.format(fget.__name__)

    @functools.wraps(fget)
    def fget_memoized(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fget(self))
        return getattr(self, attr_name)
    return property(fget_memoized)


def deprecated(message=None):
    """Decorator to mark functions as deprecated. A warning will be emitted when the function is used."""
    def deco(func):
        @functools.wraps(func)
        def wrapped(*args, **kwargs):
            warnings.warn(
                message or 'Call to deprecated function {}'.format(func.__name__),
                category=PubChemPyDeprecationWarning,
                stacklevel=2
            )
            return func(*args, **kwargs)
        return wrapped
    return deco