class Node:
    def __init__(self, func, args = None):
        self.args = [] if args is None else args
        self.func = func
    
    def get_expression(self):
        if len(self.args)==0:
            return self.func
        else:
            return self.func(*map(lambda a : a.get_expression(), self.args))

    def __str__(self):
        if len(self.args) == 0:
            return str(self.func)
        else :
            children = list(map(lambda n : str(n), self.args))
            return f'[{self.func} : {children}]'
        
    def insert(self, expr):
        if len(expr.args) == 0:
            self.args.append(
                Node(expr)
            )
        else:
            self.args.append(
                Node(expr.func)
            )
        return self.args[-1]