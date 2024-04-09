import logging
from pathlib import Path

from docutils import nodes
from sphinx.application import Sphinx as SphinxApp


_logger = logging.getLogger(f"sphinx.ext.{Path(__file__).stem}")


class InlineCE:
    def __call__(
        self,
        name: str,
        rawtext: str,
        text: str,
        lineno: int,
        options: dict = None,
        content: list = None,
    ):
        math = "".join(
            [
                r"\require{mhchem}",
                "\n",
                "\ce{" + text + "}",
            ]
        )
        node = nodes.math(rawtext, math)
        return [node], []


def setup(app: SphinxApp):

    app.add_role("ce", InlineCE())
