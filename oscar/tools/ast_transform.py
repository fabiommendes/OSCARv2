from astroid import parse
import astroid


def format_to_fstring(node):
    print(node)
    f_string_node = astroid.JoinedStr(
        lineno=node.lineno,
        col_offset=node.col_offset,
        parent=node.parent,
    )
    formatted_value_node = astroid.FormattedValue(
        lineno=node.lineno,
        col_offset=node.col_offset,
        parent=node.parent,
    )
    formatted_value_node.postinit(value=node.args[0])

    # Removes the {} since it will be represented as
    # formatted_value_node
    string = astroid.Const(node.func.expr.value.replace('{}', ''))

    f_string_node.postinit(values=[string, formatted_value_node])
    return f_string_node


def concat_to_fstring(node):
    f_string_node = astroid.JoinedStr(
        lineno=node.lineno,
        col_offset=node.col_offset,
        parent=node.parent,
    )
    formatted_value_node = astroid.FormattedValue(
        lineno=node.lineno,
        col_offset=node.col_offset,
        parent=node.parent,
    )
    formatted_value_node.postinit(value=node.args[0])

    # Removes the {} since it will be represented as
    # formatted_value_node
    string = astroid.Const(node.func.expr.value.replace('{}', ''))

    f_string_node.postinit(values=[string, formatted_value_node])
    return f_string_node


def transform(st):
    tree = parse(st)
    print(tree.as_string())


astroid.MANAGER.register_transform(astroid.Call, format_to_fstring)
astroid.MANAGER.register_transform(astroid.BinOp, concat_to_fstring)
