import astroid
from astroid import nodes
from typing import TYPE_CHECKING

from pylint.checkers import BaseChecker, IAstroidChecker

if TYPE_CHECKING:
    from pylint.lint import PyLinter


class SilentExceptionChecker(BaseChecker):
    __implements__ = IAstroidChecker

    name = "silent-exception"
    msgs = {
        "W9006": (
            "Exception without handling, only pass, detected",
            "silent-exception-handling",
            "Exception blocks should not have a bare pass",
        )
    }

    def visit_try(self, node):
        """This method is called by pylint for every try statement."""
        self._visit_try_except(node)

    def _visit_try_except(self, node):
        for handler in node.handlers:
            if self._is_silent_exception(handler):
                self.add_message("silent-exception-handling", node=handler)

    def _is_silent_exception(self, handler):
        # checks bare except
        if not handler.type:  # bare except:
            return len(handler.body) == 1 and isinstance(handler.body[0], astroid.Pass)

        # checks only pass after except
        if (
            isinstance(handler.type, astroid.Name)
            and handler.type.name in ("Exception", "BaseException")
            and len(handler.body) == 1
            and isinstance(handler.body[0], astroid.Pass)
        ):
            return True
        return False


def register(linter: "PyLinter") -> None:
    """This required method auto registers the checker during initialization."""
    linter.register_checker(SilentExceptionChecker(linter))
