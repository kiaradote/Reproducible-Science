"""Execute a notebook's code cells in order in one namespace, like Run All (stops at the first error).
No Jupyter needed. Returns a list of dicts: cell index, ok, stdout, error."""
import ast, contextlib, io, json, os, traceback

def run_notebook(path, env=None, after_cell=None, only=None):
    nb = json.load(open(path, encoding="utf-8"))
    old = {k: os.environ.get(k) for k in (env or {})}
    os.environ.update(env or {})
    ns = {"__name__": "__main__"}
    log = []
    try:
        for i, c in enumerate(nb["cells"]):
            if c["cell_type"] != "code" or (only is not None and i not in only):
                continue
            src = "".join(c["source"])
            buf = io.StringIO()
            err = None
            try:
                with contextlib.redirect_stdout(buf):
                    tree = ast.parse(src)
                    last = tree.body.pop() if tree.body and isinstance(tree.body[-1], ast.Expr) else None
                    exec(compile(tree, f"cell{i}", "exec"), ns)
                    if last is not None:
                        val = eval(compile(ast.Expression(last.value), f"cell{i}", "eval"), ns)
                        if val is not None:
                            print(val)
            except BaseException as e:
                err = f"{type(e).__name__}: {e}"
            log.append({"cell": i, "ok": err is None, "stdout": buf.getvalue(), "error": err})
            if after_cell:
                after_cell(i, ns)
            if err:
                break
    finally:
        for k, v in old.items():
            if v is None: os.environ.pop(k, None)
            else: os.environ[k] = v
    return log, ns
