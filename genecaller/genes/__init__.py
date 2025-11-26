"""Gene-specific subclasses with lazy loading."""

import importlib
from genecaller.gene import Gene

_GENE_CLASS_CACHE = {}


def get_gene_class(gene_name: str):
    """
    Get gene-specific class or fallback to base Gene class.

    Lazy loads from genecaller/genes/{gene_name.lower()}.py
    Returns cached class on subsequent calls.
    """
    if gene_name in _GENE_CLASS_CACHE:
        return _GENE_CLASS_CACHE[gene_name]

    try:
        module = importlib.import_module(f'.{gene_name.lower()}', package='genecaller.genes')
        gene_class = getattr(module, gene_name)

        if not issubclass(gene_class, Gene):
            raise TypeError(f"{gene_name} is not a subclass of Gene")

        _GENE_CLASS_CACHE[gene_name] = gene_class
        return gene_class

    except (ImportError, AttributeError, TypeError):
        _GENE_CLASS_CACHE[gene_name] = Gene
        return Gene


__all__ = ['get_gene_class', 'Gene']
