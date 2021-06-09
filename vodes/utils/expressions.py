from vodes.utils.node import Node

# Responsible for injecting noise variables
def pre_walk(node : Node, func, index=1):
    func(node, index)

    for i in range(0,len(node.args)):
        pre_walk(node.args[i], func, index * 10 + i)

def after_walk(node : Node, func, index=1):
    for i in range(0,len(node.args)):
        after_walk(node.args[i], func, index * 10 + i)
    
    func(node, index)

##
# Given a symbolic sympy expression, this function converts it into a binary tree.
##     
def convert_to_binary(expr):
    tree = convert_to_tree(expr)
    pre_walk(tree, lambda n,i : __split(n, i))
    return tree

def __split(node, index):
    if len(node.args) <= 2:
        return node

    children = [
        node.args[0],
        node.args[1]
    ]

    for i in range (2, len(node.args)):
        children[0] = Node(node.func, [children[0],children[1]])
        children[1] = node.args[i]

    node.args = children    

###
# Given a symbolic sympy expression, this function converts it into a internal tree expression.
###
def convert_to_tree(expr, tree=None):
    if tree is None:
        tree = Node(expr.func)
    
    for arg in expr.args:
        convert_to_tree(arg,tree=tree.insert(arg))

    return tree

