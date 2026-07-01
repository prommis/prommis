"""Optional Pylint plugins for PrOMMiS."""

try:
    from custom_pylint.silent_exception_checker import register as _register
except ModuleNotFoundError as exc:
    if exc.name != "custom_pylint":
        raise

    def register(_linter):
        """Fallback no-op when custom_pylint is not installed."""
        return None

else:

    def register(linter):
        """Delegate to the custom_pylint plugin when available."""
        return _register(linter)
