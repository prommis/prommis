import sys
from typing import TYPE_CHECKING
import astroid
from astroid import nodes

from pylint.checkers import BaseChecker

if TYPE_CHECKING:
    from pylint.lint import PyLinter


class SilentExceptionChecker(BaseChecker):
    name = "silent-exception"
    msgs = {
        "W9006": (
            "Exception without handling, only pass, detected. If checking if your object is a Reference, use Pyomo's built in is_reference method instead.",
            "silent-exception-handling",
            "Exception blocks should not have a bare pass",
        )
    }

    def visit_module(self, node):
        print(
            f"DEBUG: visit_module called for {node.name}", file=sys.stderr, flush=True
        )

    def visit_functiondef(self, node):
        print(
            f"DEBUG: visit_functiondef called for {node.name}",
            file=sys.stderr,
            flush=True,
        )

    def visit_tryexcept(self, node):
        import sys

        print(
            f"DEBUG: visit_tryexcept called at line {node.lineno}",
            file=sys.stderr,
            flush=True,
        )
        for handler in node.handlers:
            if self._is_silent_exception(handler):
                self.add_message("silent-exception-handling", node=handler)

    def _is_silent_exception(self, handler: astroid.ExceptHandler) -> bool:
        """Check if exception handler only contains 'pass'."""
        # Check for bare except: pass
        if not handler.type:  # bare except:
            return len(handler.body) == 1 and isinstance(handler.body[0], astroid.Pass)

        # Check for except Exception: pass or except BaseException: pass
        if (
            isinstance(handler.type, nodes.Name)
            and handler.type.name in ("Exception", "BaseException")
            and len(handler.body) == 1
            and isinstance(handler.body[0], nodes.Pass)
        ):
            return True
        return False


def register(linter):
    """Register the checker."""
    print(f"Registering plugin: {__name__}", file=sys.stderr)
    checker = SilentExceptionChecker(linter)
    linter.register_checker(checker)
    print(f"Registered checker: {checker.name}", file=sys.stderr)
