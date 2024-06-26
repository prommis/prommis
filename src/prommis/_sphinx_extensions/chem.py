import logging
from pathlib import Path
from typing import List

from docutils import nodes
from sphinx.application import Sphinx as SphinxApp
from sphinx.directives.patches import MathDirective


_logger = logging.getLogger(f"sphinx.ext.{Path(__file__).stem}")


class DirectiveCE(MathDirective):
    @property
    def equations(self) -> List[str]:
        return list(self.arguments)

    def run(self):
        lines = [
            r"\require{mhchem}",
        ]

        for eq in self.equations:
            _logger.debug("Adding equation %r", eq)
            lines.append("\ce{" + eq + "}")

        # reimplementing rest of the code from MathDirective.run()
        latex = "\n".join(lines)
        label = self.options.get("label", self.options.get("name"))
        node = nodes.math_block(
            latex,
            latex,
            classes=self.options.get("class", []),
            docname=self.env.docname,
            number=None,
            label=label,
            nowrap="nowrap" in self.options,
        )
        self.add_name(node)
        self.set_source_info(node)

        ret = [node]
        self.add_target(ret)
        return ret


class RoleCE:
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
    app.add_role("ce", RoleCE())
    app.add_directive("ce", DirectiveCE)
