# tlogo

Sequence logo generator for aligned FASTA files

## Install

```bash
pip install tlogo
```

## Development

```bash
git clone https://github.com/tmsincomb/tlogo.git
cd tlogo
pip install -e ".[dev]"
```

## Usage

```bash
tlogo --help
tlogo hello
tlogo hello YourName
```

## Testing

```bash
pytest -v
```

## Releasing

Version is managed automatically by [setuptools-scm](https://github.com/pypa/setuptools-scm) from git tags:

```bash
git tag v0.1.0
git push origin v0.1.0
```

PyPI publishing is handled by GitHub Actions via trusted publishing (OIDC). Configure your PyPI project to trust the `publish.yml` workflow.

## License

MIT
