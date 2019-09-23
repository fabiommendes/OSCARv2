import importlib as imp


def mk_all(config.mod_name):
    mod = imp.import_module(config.mod_name)
    names = [repr(x) for x in dir(mod) if not x.startswith('_')]
    names = ', '.join(names)
    return f'__all__ = [{names}]'


def mk_imports(config.mod_name):
    mod = imp.import_module(config.mod_name)
    names = [x for x in dir(mod) if not x.startswith('_')]
    return ', '.join(names)
