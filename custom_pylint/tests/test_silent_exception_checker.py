

import pytest
# from silent_exception_checker import SilentExceptionChecker
from pylint.testutils import CheckerTestCase, MessageTest
from custom_pylint.silent_exception_checker import SilentExceptionChecker

""" checks to be done: p
1a. expect: pass -> yes
1b. expect: some action -> no

2a. expect Exception: pass -> yes
2b. expect Exception: some action -> no

3a. expect BaseException: pass -> yes 
3b. expect BaseException: some action -> no 

4a. expect SpecificExcption: pass -> no
4b. expect SpecificExcption: some action -> no

5. expect Exception: pass, then some action -> no

6. no try expect -> not bad"""


class TestSilentExceptionChecker(CheckerTestCase):
    CHECKER_CLASS = SilentExceptionChecker

    def test_bare_except_pass(self):
        code = """
        try:
            foo()
        except:
            pass
        """
        node = self.get_module_node(code)
        with self.assertAddsMessages(
            MessageTest(
                msg_id="silent-exception-handling",
                node=node.body[0].handlers[0],
            )
        ):
            self.walk(node)

    def test_bare_except_exception_with_code(self):
        code = """
        try:
            foo()
        except Exception:
            print("error")
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)

    def test_except_exception_pass(self):
        code = """
        try:
            foo()
        except Exception:
            pass
        """
        node = self.get_module_node(code)
        with self.assertAddsMessages(
            MessageTest(
                msg_id="silent-exception-handling",
                node=node.body[0].handlers[0],
            )
        ):
            self.walk(node)

    def test_except_exception_with_code(self):
        code = """
        try:
            foo()
        except Exception:
            print("error")
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)

    def test_except_baseexception_pass(self):
        code = """
        try:
            foo()
        except BaseException:
            pass
        """
        node = self.get_module_node(code)
        with self.assertAddsMessages(
            MessageTest(
                msg_id="silent-exception-handling",
                node=node.body[0].handlers[0],
            )
        ):
            self.walk(node)


    def test_except_baseexception_with_code(self):
        code = """
        try:
            foo()
        except Exception:
            print("error")
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)

    def test_except_specific_exception_pass(self):
        code = """
        try:
            foo()
        except ValueError:
            pass
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)

    def test_except_specific_exception_with_code(self):
        code = """
        try:
            foo()
        except ValueError:
            print("error")
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)


    def test_except_with_multiple_statements(self):
        code = """
        try:
            foo()
        except Exception:
            pass
            print("error")
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)

    def test_no_tryexcept(self):
        code = """
            foo()
        """
        node = self.get_module_node(code)
        with self.assertNoMessages():
            self.walk(node)