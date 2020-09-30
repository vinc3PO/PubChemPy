import json


class PubChemPyDeprecationWarning(Warning):
    """Warning category for deprecated features."""
    pass


class PubChemPyError(Exception):
    """Base class for all PubChemPy exceptions."""
    pass


class ResponseParseError(PubChemPyError):
    """PubChem response is uninterpretable."""
    pass


class PubChemHTTPError(PubChemPyError):
    """Generic error class to handle all HTTP error codes."""
    def __init__(self, e):
        self.code = e.code
        self.msg = e.reason
        try:
            self.msg += ': %s' % json.loads(e.read().decode())['Fault']['Details'][0]
        except (ValueError, IndexError, KeyError):
            pass
        if self.code == 400:
            raise BadRequestError(self.msg)
        elif self.code == 404:
            raise NotFoundError(self.msg)
        elif self.code == 405:
            raise MethodNotAllowedError(self.msg)
        elif self.code == 504:
            raise TimeoutError(self.msg)
        elif self.code == 501:
            raise UnimplementedError(self.msg)
        elif self.code == 500:
            raise ServerError(self.msg)

    def __str__(self):
        return repr(self.msg)


class BadRequestError(PubChemHTTPError):
    """Request is improperly formed (syntax error in the URL, POST body, etc.)."""
    def __init__(self, msg='Request is improperly formed'):
        self.msg = msg


class NotFoundError(PubChemHTTPError):
    """The input record was not found (e.g. invalid CID)."""
    def __init__(self, msg='The input record was not found'):
        self.msg = msg


class MethodNotAllowedError(PubChemHTTPError):
    """Request not allowed (such as invalid MIME type in the HTTP Accept header)."""
    def __init__(self, msg='Request not allowed'):
        self.msg = msg


class TimeoutError(PubChemHTTPError):
    """The request timed out, from server overload or too broad a request.

    See :ref:`Avoiding TimeoutError <avoiding_timeouterror>` for more information.
    """
    def __init__(self, msg='The request timed out'):
        self.msg = msg


class UnimplementedError(PubChemHTTPError):
    """The requested operation has not (yet) been implemented by the server."""
    def __init__(self, msg='The requested operation has not been implemented'):
        self.msg = msg


class ServerError(PubChemHTTPError):
    """Some problem on the server side (such as a database server down, etc.)."""
    def __init__(self, msg='Some problem on the server side'):
        self.msg = msg
