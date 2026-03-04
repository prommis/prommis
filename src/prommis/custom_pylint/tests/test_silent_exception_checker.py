import pytest
import astroid
import pylint.testutils
from prommis.custom_pylint.silent_exception_checker import SilentExceptionChecker


class TestSilentExceptionChecker(pylint.testutils.CheckerTestCase):
    CHECKER_CLASS = SilentExceptionChecker

    # ------------------------------------------------------------------
    # POSITIVE tests — checker SHOULD emit a message
    # ------------------------------------------------------------------

    def test_bare_except_with_pass_is_flagged(self):
        """bare except: pass must be flagged."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except:  #@
                pass
        """)
        with self.assertAddsMessages(
            pylint.testutils.MessageTest(
                msg_id="silent-exception-handling",
                node=handler_node,
            )
        ):
            self.checker.visit_tryexcept(handler_node.parent)

    def test_except_exception_with_pass_is_flagged(self):
        """except Exception: pass must be flagged."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except Exception:  #@
                pass
        """)
        with self.assertAddsMessages(
            pylint.testutils.MessageTest(
                msg_id="silent-exception-handling",
                node=handler_node,
            )
        ):
            self.checker.visit_tryexcept(handler_node.parent)

    def test_except_base_exception_with_pass_is_flagged(self):
        """except BaseException: pass must be flagged."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except BaseException:  #@
                pass
        """)
        with self.assertAddsMessages(
            pylint.testutils.MessageTest(
                msg_id="silent-exception-handling",
                node=handler_node,
            )
        ):
            self.checker.visit_tryexcept(handler_node.parent)

    def test_except_exception_with_alias_and_pass_is_flagged(self):
        """except Exception as e: pass must also be flagged."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except Exception as e:  #@
                pass
        """)
        with self.assertAddsMessages(
            pylint.testutils.MessageTest(
                msg_id="silent-exception-handling",
                node=handler_node,
            )
        ):
            self.checker.visit_tryexcept(handler_node.parent)

    # ------------------------------------------------------------------
    # NEGATIVE tests — checker should NOT emit any message
    # ------------------------------------------------------------------

    def test_except_exception_with_logging_not_flagged(self):
        """except Exception with a log call is not silent — no message."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except Exception:  #@
                print("error")
        """)
        with self.assertNoMessages():
            self.checker.visit_tryexcept(handler_node.parent)

    def test_except_exception_with_raise_not_flagged(self):
        """except Exception: raise is fine — not a silent swallow."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except Exception:  #@
                raise
        """)
        with self.assertNoMessages():
            self.checker.visit_tryexcept(handler_node.parent)

    def test_except_specific_exception_with_pass_not_flagged(self):
        """except ValueError: pass should NOT be flagged (checker only targets broad exceptions)."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except ValueError:  #@
                pass
        """)
        with self.assertNoMessages():
            self.checker.visit_tryexcept(handler_node.parent)

    def test_except_exception_with_multiple_statements_not_flagged(self):
        """More than one statement in the body — not flagged."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except Exception:  #@
                pass
                print("done")
        """)
        with self.assertNoMessages():
            self.checker.visit_tryexcept(handler_node.parent)

    def test_bare_except_with_non_pass_body_not_flagged(self):
        """Bare except with a real body — not flagged."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except:  #@
                x = 0
        """)
        with self.assertNoMessages():
            self.checker.visit_tryexcept(handler_node.parent)

    # ------------------------------------------------------------------
    # KNOWN GAPS — these tests document current checker LIMITATIONS.
    # They are marked xfail so CI stays green while tracking the gaps.
    # Remove xfail once the checker is extended to cover these cases.
    # ------------------------------------------------------------------

    @pytest.mark.xfail(
        reason="Checker does not yet detect tuple-typed except clauses "
               "that include Exception, e.g. except (ValueError, Exception): pass"
    )
    def test_except_tuple_including_exception_with_pass_is_flagged(self):
        """except (ValueError, Exception): pass — broad catch hidden in a tuple."""
        (handler_node,) = astroid.extract_node("""
        def f():
            try:
                x = 1
            except (ValueError, Exception):  #@
                pass
        """)
        with self.assertAddsMessages(
            pylint.testutils.MessageTest(
                msg_id="silent-exception-handling",
                node=handler_node,
            )
        ):
            self.checker.visit_tryexcept(handler_node.parent)

    @pytest.mark.xfail(
        reason="Checker does not yet handle qualified names like builtins.Exception"
    )
    def test_except_qualified_exception_with_pass_is_flagged(self):
        """except builtins.Exception: pass — qualified broad catch."""
        (handler_node,) = astroid.extract_node("""
        import builtins
        def f():
            try:
                x = 1
            except builtins.Exception:  #@
                pass
        """)
        with self.assertAddsMessages(
            pylint.testutils.MessageTest(
                msg_id="silent-exception-handling",
                node=handler_node,
            )
        ):
            self.checker.visit_tryexcept(handler_node.parent)